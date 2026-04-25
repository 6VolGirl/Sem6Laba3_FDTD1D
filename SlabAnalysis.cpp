//
// Created by 6anna on 25.04.2026.
//


#include "SlabAnalysis.h"


#include <fstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>




// inc(t), tot(t), eref(t), trans(t)
static void writeSlabTimeSeriesCSV(const std::string&   filename,
                                    const FieldMonitor&  incMon,
                                    const FieldMonitor&  totMon,
                                    const FieldMonitor&  transMon,
                                    double               dt)
{
    const int N = static_cast<int>(
    std::min({ incMon.dataEx.size(),
               incMon.dataHy.size(),
               totMon.dataEx.size(),
               totMon.dataHy.size(),
               transMon.dataEx.size(),
               transMon.dataHy.size() })
    );

    if (N == 0) {
        std::cerr << "writeSlabTimeSeriesCSV: нет данных, файл " << filename << " не записан\n";
        return;
    }

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "writeSlabTimeSeriesCSV: не удалось открыть " << filename << "\n";
        return;
    }

    out << "step,time,"
        << "Ex_inc,Hy_inc,S_inc,"
        << "Ex_tot,Hy_tot,S_tot,"
        << "Ex_ref,Hy_ref,S_ref,"
        << "Ex_trans,Hy_trans,S_trans\n";
    out << std::scientific << std::setprecision(10);

    for (int n = 0; n < N; ++n) {
        const double t    = n * dt;

        const double exInc = incMon.dataEx[n];
        const double hyInc = incMon.dataHy[n];
        const double sInc  = exInc * hyInc;

        const double exTot = totMon.dataEx[n];
        const double hyTot = totMon.dataHy[n];
        const double sTot  = exTot * hyTot;

        const double exRef = exTot - exInc;
        const double hyRef = hyTot - hyInc;
        const double sRef  = exRef * hyRef;

        const double exTr  = transMon.dataEx[n];
        const double hyTr  = transMon.dataHy[n];
        const double sTr   = exTr * hyTr;

        out << n     << "," << t     << ","
            << exInc << "," << hyInc << "," << sInc << ","
            << exTot << "," << hyTot << "," << sTot << ","
            << exRef << "," << hyRef << "," << sRef << ","
            << exTr  << "," << hyTr  << "," << sTr  << "\n";
    }

    std::cout << "Timeseries written to " << filename
              << " (" << N << " steps)\n";
}


std::vector<SlabResult> SlabAnalysis::runSlab(double L_nm, int numSteps)
{
    const int slabCenter = base_.nx / 2;
    const int halfL      = static_cast<int>(std::round(L_nm / cfgS_.a_nm / 2.0));
    const int slabStart  = slabCenter - halfL;
    const int slabEnd    = slabCenter + halfL;

    if (slabStart <= base_.source_pos)
        throw std::runtime_error("SlabAnalysis: пластина перекрывает источник");
    if (slabEnd + cfgS_.monTransOffset >= base_.nx - base_.pmlThickness)
        throw std::runtime_error("SlabAnalysis: монитор T слишком близко к PML");

    const int monRefPos   = base_.source_pos - cfgS_.monRefOffset;
    const int monTransPos = slabEnd + cfgS_.monTransOffset;

    if (monRefPos <= base_.pmlThickness)
        throw std::runtime_error("SlabAnalysis: монитор R слишком близко к PML");

    std::cout << "SlabAnalysis L=" << L_nm << " nm:\n"
              << "  slabStart=" << slabStart
              << "  slabEnd="   << slabEnd   << "\n"
              << "  monRef="    << monRefPos
              << "  monTrans="  << monTransPos << "\n";

    SimulationParameters pRun = base_;

    if (numSteps > 0)
        pRun.numTimeSteps = numSteps;

    Monitor      mon(pRun.dt, monRefPos);
    FieldMonitor transMon(monTransPos);


    {
        FDTD1D sim(pRun);
        sim.attachMonitor(&mon.incMonitor());
        sim.attachMonitor(&transMon);
        sim.run();
    }

    const std::vector<double> incRefEx   = mon.incMonitor().dataEx;
    const std::vector<double> incRefHy   = mon.incMonitor().dataHy;
    const std::vector<double> transIncEx = transMon.dataEx;
    const std::vector<double> transIncHy = transMon.dataHy;

    transMon.clear();

     {
        MaterialLayout layout;
        layout.add(MaterialRegion::Silica(slabStart, slabEnd));

        FDTD1D simWork(pRun, layout);
        simWork.attachMonitor(&mon.totMonitor());
        simWork.attachMonitor(&transMon);
        simWork.run();

        simWork.writeFieldCSV("slab_L" + std::to_string(static_cast<int>(L_nm)) + "_field.csv");
    }

    writeSlabTimeSeriesCSV(
        "slab_L" + std::to_string(static_cast<int>(L_nm)) + "nm_monitor.csv",
        mon.incMonitor(),
        mon.totMonitor(),
        transMon,
        pRun.dt
    );

    const double fMin = cfgS_.a_nm / cfgS_.lambdaMax;
    const double fMax = cfgS_.a_nm / cfgS_.lambdaMin;
    const int    nF   = cfgS_.nWavelengths;


    const int Nref = static_cast<int>(
        std::min(incRefEx.size(), mon.totMonitor().dataEx.size())
    );

    const int Ntr = static_cast<int>(
        std::min(transIncEx.size(), transMon.dataEx.size())
    );

    if (Nref == 0 || Ntr == 0)
        throw std::runtime_error("SlabAnalysis: недостаточно данных на мониторах");

    std::vector<double> eInc(Nref), eRef(Nref);
    for (int n = 0; n < Nref; ++n) {
        eInc[n] = incRefEx[n];
        eRef[n] = mon.totMonitor().dataEx[n] - incRefEx[n];
    }

    std::vector<double> eTransInc(Ntr), eTransWork(Ntr);
    for (int n = 0; n < Ntr; ++n) {
        eTransInc[n]  = transIncEx[n];
        eTransWork[n] = transMon.dataEx[n];
    }

    auto dftManual = [&](const std::vector<double>& sig, int N, double freq) -> double
    {
        double re = 0.0, im = 0.0;
        const double omega = 2.0 * M_PI * freq;
        for (int n = 0; n < N; ++n) {
            const double phase = -omega * n * pRun.dt;
            re += sig[n] * std::cos(phase);
            im += sig[n] * std::sin(phase);
        }
        return std::sqrt(re * re + im * im);
    };

    std::vector<double> freqs(nF);
    const double df = (nF > 1) ? (fMax - fMin) / (nF - 1) : 0.0;
    for (int i = 0; i < nF; ++i)
        freqs[i] = fMin + i * df;

    double maxAmpIncRef = 0.0;
    double maxAmpIncTr  = 0.0;
    for (int i = 0; i < nF; ++i) {
        maxAmpIncRef = std::max(maxAmpIncRef, dftManual(eInc,      Nref, freqs[i]));
        maxAmpIncTr  = std::max(maxAmpIncTr,  dftManual(eTransInc, Ntr,  freqs[i]));
    }

    const double thresholdRef = maxAmpIncRef * 1e-3;
    const double thresholdTr  = maxAmpIncTr  * 1e-3;

    result_.clear();
    result_.reserve(nF);

    for (int i = nF - 1; i >= 0; --i) {
        const double f   = freqs[i];
        const double lam = (f > 1e-30) ? (cfgS_.a_nm / f) : 0.0;

        // R = |E_ref|^2 / |E_inc|^2
        const double ampIncRef = dftManual(eInc, Nref, f);
        const double ampRef    = dftManual(eRef, Nref, f);

        double R = 0.0;
        if (ampIncRef > thresholdRef) {
            R = (ampRef * ampRef) / (ampIncRef * ampIncRef);
        }


        // T = |E_trans|^2 / |E_inc(trans)|^2
        const double ampIncTr   = dftManual(eTransInc,  Ntr, f);
        const double ampTrans   = dftManual(eTransWork, Ntr, f);

        double T = 0.0;
        if (ampIncTr > thresholdTr) {
            T = (ampTrans * ampTrans) / (ampIncTr * ampIncTr);
        }

        const double R_th = theoreticalR(lam, L_nm, cfgS_.n1, cfgS_.n2);
        const double T_th = theoreticalT(lam, L_nm, cfgS_.n1, cfgS_.n2);

        result_.push_back({lam, R, T, R_th, T_th});
    }

    return result_;
}


void SlabAnalysis::runMultipleL(const std::vector<double>& L_vals,
                                const std::string& prefix)
{
    for (double L : L_vals) {
        std::cout << "\n=== Slab L=" << L << " nm ===\n";
        runSlab(L);
        writeCSV(prefix + std::to_string(static_cast<int>(L)) + "nm.csv",
                 result_);
    }
}

double SlabAnalysis::theoreticalR(double lambda_nm, double L_nm,
                                   double n1, double n2)
{
    // Формула Фабри-Перо для пластины толщиной L_nm (нормальное падение)
    const double r = (n1 - n2) / (n1 + n2);  // амплитудный коэф. отражения на границе
    const double r2 = r * r;
    // Фаза за двойной проход внутри пластины
    const double phi = 2.0 * M_PI * n2 * L_nm / lambda_nm;  // lambda в тех же единицах
    const double sinPhi2 = std::sin(phi) * std::sin(phi);
    // R = r^2 * 4sin^2(phi) / (1 - r^4 + 4r^2 sin^2(phi))  — полная формула
    // Упрощение Айри:
    const double num   = 4.0 * r2 * sinPhi2;
    const double denom = (1.0 - r2) * (1.0 - r2) + 4.0 * r2 * sinPhi2;
    return (denom > 1e-30) ? num / denom : 0.0;
}

double SlabAnalysis::theoreticalT(double lambda_nm, double L_nm,
                                   double n1, double n2)
{
    return 1.0 - theoreticalR(lambda_nm, L_nm, n1, n2);
}

void SlabAnalysis::writeCSV(const std::string& filename,
                             const std::vector<SlabResult>& data)
{
    if (data.empty()) {
        std::cerr << "SlabAnalysis::writeCSV: нет данных\n";
        return;
    }
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "SlabAnalysis::writeCSV: не удалось открыть " << filename << "\n";
        return;
    }
    out << "lambda_nm,R_fdtd,T_fdtd,R_theory,T_theory\n";
    out << std::scientific << std::setprecision(10);
    for (const auto& row : data)
        out << row.lambda_nm << ","
            << row.R_fdtd   << ","
            << row.T_fdtd   << ","
            << row.R_theory << ","
            << row.T_theory << "\n";

    std::cout << "SlabAnalysis: written " << data.size()
              << " points to " << filename << "\n";
}