//
// Created by 6anna on 26.04.2026.
//

#include "PhotonCrystalAnalysis.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iomanip>


PhotonicCrystalAnalysis::PhotonicCrystalAnalysis(const SimulationParameters& base)
    : base_(base), cfg_() {}

PhotonicCrystalAnalysis::PhotonicCrystalAnalysis(const SimulationParameters& base,
                                                 const ConfigPhot& cfg)
    : base_(base), cfg_(cfg) {}

double PhotonicCrystalAnalysis::dftAmp(const std::vector<double>& sig,
                                        int N, double freq, double dt)
{
    double re = 0.0, im = 0.0;
    const double omega = 2.0 * M_PI * freq;
    for (int n = 0; n < N; ++n) {
        const double phase = -omega * n * dt;
        re += sig[n] * std::cos(phase);
        im += sig[n] * std::sin(phase);
    }
    return std::sqrt(re * re + im * im);
}

MaterialLayout PhotonicCrystalAnalysis::buildLayout(int crystalStart,
                                            int dA_cells, int dB_cells,
                                                 int numPeriods)
{
    MaterialLayout layout;
    int pos = crystalStart;
    for (int p = 0; p < numPeriods; ++p) {
        layout.add(MaterialRegion::TiO2(pos, pos + dA_cells));
        pos += dA_cells;

        layout.add(MaterialRegion::Silica(pos, pos + dB_cells));
        pos += dB_cells;
    }
    return layout;
}


std::vector<PCResult> PhotonicCrystalAnalysis::runOnce(int crystalStart,
                                                       int dA_cells, int dB_cells,
                                                       int numPeriods) const
{
    const int crystalEnd = crystalStart + numPeriods * (dA_cells + dB_cells);

    const int monRefPos   = base_.source_pos - cfg_.monRefOffset;
    const int monTransPos = crystalEnd + cfg_.monTransOffset;

    if (monRefPos <= base_.pmlThickness)
        throw std::runtime_error("PC: монитор R слишком близко к PML");
    if (monTransPos >= base_.nx - base_.pmlThickness)
        throw std::runtime_error("PC: монитор T слишком близко к PML");
    if (crystalStart <= base_.source_pos)
        throw std::runtime_error("PC: кристалл перекрывает источник");

    std::cout << "  crystalStart=" << crystalStart
              << "  crystalEnd="   << crystalEnd
              << "  monRef="       << monRefPos
              << "  monTrans="     << monTransPos << "\n";

    SimulationParameters pRun = base_;

    Monitor      mon(pRun.dt, monRefPos);
    FieldMonitor transMon(monTransPos);

    // вакуум
    {
        FDTD1D sim(pRun);
        sim.attachMonitor(&mon.incMonitor());
        sim.attachMonitor(&transMon);
        sim.run();
    }

    const std::vector<double> eIncRef(mon.incMonitor().dataEx);
    const std::vector<double> eTransInc(transMon.dataEx);
    transMon.clear();

    // с кристаллом
    {
        MaterialLayout layout = buildLayout(crystalStart, dA_cells,
                                            dB_cells, numPeriods);
        FDTD1D simWork(pRun, layout);
        simWork.attachMonitor(&mon.totMonitor());
        simWork.attachMonitor(&transMon);
        simWork.run();
    }

    const double fMin = cfg_.a_nm / cfg_.lambdaMax;
    const double fMax = cfg_.a_nm / cfg_.lambdaMin;
    const int    nF   = cfg_.nWavelengths;

    const int Nref = static_cast<int>(std::min(eIncRef.size(), mon.totMonitor().dataEx.size()));
    const int Ntr  = static_cast<int>( std::min(eTransInc.size(), transMon.dataEx.size()));

    // Построить E_ref = E_tot - E_inc
    std::vector<double> eRef(Nref), eInc(Nref);
    for (int n = 0; n < Nref; ++n) {
        eInc[n] = eIncRef[n];
        eRef[n] = mon.totMonitor().dataEx[n] - eIncRef[n];
    }

    std::vector<double> eTrInc(Ntr), eTrWork(Ntr);
    for (int n = 0; n < Ntr; ++n) {
        eTrInc[n]  = eTransInc[n];
        eTrWork[n] = transMon.dataEx[n];
    }

    // Найти максимумы для порогов
    double maxIncRef = 0.0, maxIncTr = 0.0;
    const double df = (nF > 1) ? (fMax - fMin) / (nF - 1) : 0.0;
    for (int i = 0; i < nF; ++i) {
        double f = fMin + i * df;
        maxIncRef = std::max(maxIncRef, dftAmp(eInc,   Nref, f, pRun.dt));
        maxIncTr  = std::max(maxIncTr,  dftAmp(eTrInc, Ntr,  f, pRun.dt));
    }
    const double thrRef = maxIncRef * 1e-3;
    const double thrTr  = maxIncTr  * 1e-3;

    std::vector<PCResult> result;
    result.reserve(nF);

    for (int i = nF - 1; i >= 0; --i) {
        const double f   = fMin + i * df;
        const double lam = (f > 1e-30) ? (cfg_.a_nm / f) : 0.0;

        // R = |E_ref|^2 / |E_inc|^2
        const double ampInc = dftAmp(eInc,   Nref, f, pRun.dt);
        const double ampRef = dftAmp(eRef,   Nref, f, pRun.dt);
        const double R = (ampInc > thrRef) ? (ampRef * ampRef) / (ampInc * ampInc) : 0.0;

        // T = |E_trans|^2 / |E_inc_trans|^2
        const double ampTrInc  = dftAmp(eTrInc,  Ntr, f, pRun.dt);
        const double ampTrWork = dftAmp(eTrWork, Ntr, f, pRun.dt);
        const double T = (ampTrInc > thrTr) ? (ampTrWork * ampTrWork) / (ampTrInc * ampTrInc) : 0.0;

        result.push_back({lam, R, T, R + T});
    }

    return result;
}


void PhotonicCrystalAnalysis::runEqualThickness(const std::vector<int>& periodCounts,
                                                const std::string& outPrefix)
{
    const int dA = static_cast<int>(std::round(cfg_.d_A_equal / cfg_.a_nm));
    const int dB = static_cast<int>(std::round(cfg_.d_B_equal / cfg_.a_nm));

    // Начало кристалла — правее источника
    const int crystalStart = base_.source_pos + cfg_.monRefOffset + 50;

    for (int N : periodCounts) {
        std::cout << "\n=== PC equal thickness, N=" << N
                  << " (TiO2=" << cfg_.d_A_equal
                  << " nm, SiO2=" << cfg_.d_B_equal << " nm) ===\n";

        auto res = runOnce(crystalStart, dA, dB, N);

        const std::string fname = outPrefix + "equal_N" + std::to_string(N) + ".csv";
        writeCSV(fname, res);
    }
}

void PhotonicCrystalAnalysis::runEqualOptPath(const std::vector<int>& periodCounts,
                                              const std::string& outPrefix)
{
    // d_A * n_A = d_B * n_B → четвертьволновый стек
    const int dA = static_cast<int>(std::round(cfg_.d_A_opt / cfg_.a_nm));
    const int dB = static_cast<int>(std::round(cfg_.d_B_opt / cfg_.a_nm));

    const int crystalStart = base_.source_pos + cfg_.monRefOffset + 50;

    for (int N : periodCounts) {
        std::cout << "\n=== PC equal optical path, N=" << N
                  << " (TiO2=" << cfg_.d_A_opt
                  << " nm, SiO2=" << cfg_.d_B_opt << " nm) ===\n";

        auto res = runOnce(crystalStart, dA, dB, N);

        const std::string fname = outPrefix + "optpath_N" + std::to_string(N) + ".csv";
        writeCSV(fname, res);
    }
}

void PhotonicCrystalAnalysis::writeCSV(const std::string& filename,
                                        const std::vector<PCResult>& data)
{
    if (data.empty()) {
        std::cerr << "PhotonicCrystalAnalysis::writeCSV: нет данных\n";
        return;
    }
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "PhotonicCrystalAnalysis: не открыть " << filename << "\n";
        return;
    }
    out << "lambda_nm,R_fdtd,T_fdtd,RT_balance\n";
    out << std::scientific << std::setprecision(10);
    for (const auto& row : data)
        out << row.lambda_nm << ","
            << row.R_fdtd    << ","
            << row.T_fdtd    << ","
            << row.balance   << "\n";

    std::cout << "PC: written " << data.size()
              << " points to " << filename << "\n";
}
