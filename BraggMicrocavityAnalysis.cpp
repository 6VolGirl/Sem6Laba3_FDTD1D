//
// Created by 6anna on 26.04.2026.
//

#include "BraggMicrocavityAnalysis.h"



#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <complex>
#include <numeric>

// ---------------------------------------------------------------------------
// Геометрия слоёв
// ---------------------------------------------------------------------------
void CavityAnalysis::computeLayerCells(int& dH, int& dL, int& dCav) const
{
    dH = static_cast<int>(std::round(cfgC_.lambda0_nm / (4.0 * cfgC_.nH   * cfgC_.dx_nm)));
    dL = static_cast<int>(std::round(cfgC_.lambda0_nm / (4.0 * cfgC_.nL   * cfgC_.dx_nm)));

    if (cfgC_.cavThickness_nm > 0.0)
        dCav = static_cast<int>(std::round(cfgC_.cavThickness_nm / cfgC_.dx_nm));
    else
        dCav = static_cast<int>(std::round(cfgC_.lambda0_nm / (2.0 * cfgC_.nCav * cfgC_.dx_nm)));

    std::cout << "CavityAnalysis geometry:\n"
              << "  dH   = " << dH   << " cells (" << dH   * cfgC_.dx_nm << " nm)\n"
              << "  dL   = " << dL   << " cells (" << dL   * cfgC_.dx_nm << " nm)\n"
              << "  dCav = " << dCav << " cells (" << dCav * cfgC_.dx_nm << " nm)\n"
              << "  nPairs = " << cfgC_.nPairs << "\n";
}

// ---------------------------------------------------------------------------
// Вспомогательная: ДПФ потока Пойнтинга S(f) = Re[E(f)*conj(H(f))] / 2
// Принимает временные ряды Ex и Hy уже записанные FieldMonitor::record()
// ---------------------------------------------------------------------------
static std::complex<double> dftOne(const std::vector<double>& sig,
                                   double freq, double dt)
{
    std::complex<double> acc(0.0, 0.0);
    const double omega = 2.0 * M_PI * freq;
    const int N = static_cast<int>(sig.size());
    for (int n = 0; n < N; ++n) {
        const double ph = -omega * n * dt;
        acc += sig[n] * std::complex<double>(std::cos(ph), std::sin(ph));
    }
    return acc;
}

static double poyntingFlux(const std::vector<double>& ex,
                            const std::vector<double>& hy,
                            double freq, double dt)
{
    // S(f) = Re[ E_x(f) * conj(H_y(f)) ] / 2
    return 0.5 * std::real(dftOne(ex, freq, dt) * std::conj(dftOne(hy, freq, dt)));
}

// ---------------------------------------------------------------------------
// Вакуумный прогон
// FDTD1D::attachMonitor принимает один FieldMonitor*, поэтому делаем
// два прогона — один для позиции monRefPos, другой для monTransPos.
// ---------------------------------------------------------------------------
void CavityAnalysis::runVacuum(int monRefPos, int monTransPos,
                               FieldMonitor& incRef,
                               FieldMonitor& incTrans)
{
    SimulationParameters pVac = base_;
    pVac.injectionType = SimulationParameters::SOFT;

    {
        FDTD1D sim(pVac);
        sim.attachMonitor(&incRef);
        sim.run();
    }
    {
        FDTD1D sim(pVac);
        sim.attachMonitor(&incTrans);
        sim.run();
    }

    std::cout << "  Vacuum done. incRef=" << incRef.size()
              << "  incTrans=" << incTrans.size() << "\n";
}

// ---------------------------------------------------------------------------
// Один прогон со структурой
//
// R(f) через поток Пойнтинга:
//   E_ref = E_tot - E_inc   (на монitorе отражения)
//   H_ref = H_tot - H_inc
//   S_ref(f)     = Re[E_ref(f)*conj(H_ref(f))] / 2  — всегда <= 0 (волна в -x)
//   S_inc_ref(f) = Re[E_inc(f)*conj(H_inc(f))] / 2  — всегда >= 0 (волна в +x)
//   R(f) = -S_ref(f) / S_inc_ref(f)
//
// T(f) через поток Пойнтинга:
//   S_trans(f)     = Re[E_trans(f)*conj(H_trans(f))] / 2
//   S_inc_trans(f) = Re[E_inc_trans(f)*conj(H_inc_trans(f))] / 2
//   T(f) = S_trans(f) / S_inc_trans(f)
//
// Monitor::computeTransmission() уже реализован для T.
// Здесь R тоже пересчитываем через Пойнтинга, а не через computeReflection().
// ---------------------------------------------------------------------------
void CavityAnalysis::runOnce(const MaterialLayout& layout,
                              int monRefPos, int monTransPos,
                              bool saveField,
                              const FieldMonitor& incRef,
                              const FieldMonitor& incTrans,
                              std::vector<double>& outR,
                              std::vector<double>& outT,
                              std::vector<double>& outFreqs)
{
    SimulationParameters pWork = base_;
    pWork.injectionType = SimulationParameters::SOFT;

    Monitor mon(pWork.dt, monRefPos);

    // ── ИСПРАВЛЕНИЕ: копируем вакуумный прогон напрямую в incMon_,
    //    а не запускаем лишний simInc (который давал E_ref = 0).
    mon.incMonitor().dataEx = incRef.dataEx;
    mon.incMonitor().dataHy = incRef.dataHy;

    // Прогон total (со структурой) на мониторе отражения
    {
        FDTD1D simRef(pWork, layout);
        simRef.attachMonitor(&mon.totMonitor());
        simRef.run();
    }

    // Прогон transmitted (со структурой)
    FieldMonitor transWork(monTransPos);
    {
        FDTD1D simT(pWork, layout);
        simT.attachMonitor(&transWork);
        simT.run();

        if (saveField)
            simT.writeFieldCSV(cfgC_.fieldPrefix + "_raw.csv");
    }

    const double fMin = cfgC_.dx_nm / cfgC_.lambdaMax;
    const double fMax = cfgC_.dx_nm / cfgC_.lambdaMin;
    const int nF = cfgC_.nWavelengths;

    outFreqs.resize(nF);
    outR.resize(nF);
    outT.resize(nF);

    auto specR = mon.computeReflectionPoynting(fMin, fMax, nF);
    auto specT = mon.computeTransmissionPoynting(transWork, incTrans, fMin, fMax, nF);

    for (int i = 0; i < nF; ++i) {
        outFreqs[i] = specR[i].first;
        outR[i]     = specR[i].second;
        outT[i]     = specT[i].second;
    }
}

// ---------------------------------------------------------------------------
// Главная функция
// ---------------------------------------------------------------------------
std::vector<CavityResult> CavityAnalysis::runCavity()
{
    int dH, dL, dCav;
    computeLayerCells(dH, dL, dCav);

    const int mirrorCells = cfgC_.nPairs * (dH + dL);
    const int totalStruct = 2 * mirrorCells + dCav;

    const int minStart = base_.source_pos + 150;
    const int maxEnd   = base_.nx - base_.pmlThickness - cfgC_.monTransOffset - 20;

    if (maxEnd - minStart < totalStruct)
        throw std::runtime_error(
            "CavityAnalysis: структура не помещается в сетку; увеличьте nx или уменьшите nPairs");

    structStart_ = minStart + ((maxEnd - minStart) - totalStruct) / 2;
    structEnd_   = structStart_ + totalStruct;
    cavStart_    = structStart_ + mirrorCells;
    cavEnd_      = cavStart_ + dCav;

    const int monRefPos   = base_.source_pos - cfgC_.monRefOffset;
    const int monTransPos = structEnd_ + cfgC_.monTransOffset;

    if (monRefPos <= 1 || monTransPos >= base_.nx - 1)
        throw std::runtime_error("CavityAnalysis: мониторы выходят за границы сетки");

    std::cout << "CavityAnalysis::runCavity()\n"
              << "  structStart=" << structStart_
              << "  cavStart="    << cavStart_
              << "  cavEnd="      << cavEnd_
              << "  structEnd="   << structEnd_  << "\n"
              << "  monRef="      << monRefPos
              << "  monTrans="    << monTransPos << "\n";

    // Вакуумный прогон
    FieldMonitor incRef(monRefPos);
    FieldMonitor incTrans(monTransPos);
    runVacuum(monRefPos, monTransPos, incRef, incTrans);

    // Полость (DBR | cavity | DBR)
    MaterialLayout layoutCav;
    layoutCav.addBraggCavity(structStart_, cfgC_.nH, cfgC_.nL,
                             dH, dL, cfgC_.nPairs, dCav);

    std::vector<double> Rcav, Tcav, freqs;
    runOnce(layoutCav, monRefPos, monTransPos, true,
            incRef, incTrans, Rcav, Tcav, freqs);

    // Периодическая структура (HL)^{2N+1}
    MaterialLayout layoutPer;
    {
        int cur = structStart_;
        const int totalPairs = 2 * cfgC_.nPairs + 1;
        for (int p = 0; p < totalPairs; ++p) {
            layoutPer.add(MaterialRegion(cur, cur + dH, cfgC_.nH * cfgC_.nH, 1.0, "H"));
            cur += dH;
            layoutPer.add(MaterialRegion(cur, cur + dL, cfgC_.nL * cfgC_.nL, 1.0, "L"));
            cur += dL;
        }
    }

    std::vector<double> Rper, Tper, freqsPer;
    runOnce(layoutPer, monRefPos, monTransPos, false,
            incRef, incTrans, Rper, Tper, freqsPer);

    const int nF = cfgC_.nWavelengths;
    result_.clear();
    result_.reserve(nF);

    for (int i = nF - 1; i >= 0; --i) {
        const double f   = freqs[i];
        const double lam = (f > 1e-30) ? (cfgC_.dx_nm / f) : 0.0;

        CavityResult r{};
        r.lambda_nm   = lam;
        r.R_cavity    = (i < (int)Rcav.size()) ? Rcav[i] : 0.0;
        r.T_cavity    = (i < (int)Tcav.size()) ? Tcav[i] : 0.0;
        r.R_periodic  = (i < (int)Rper.size()) ? Rper[i] : 0.0;
        r.T_periodic  = (i < (int)Tper.size()) ? Tper[i] : 0.0;
        result_.push_back(r);
    }

    return result_;
}

// ---------------------------------------------------------------------------
// CSV экспорт
// ---------------------------------------------------------------------------
void CavityAnalysis::writeSpectraCSV(const std::string& filename) const
{
    if (result_.empty()) {
        std::cerr << "CavityAnalysis::writeSpectraCSV: нет данных\n";
        return;
    }

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "CavityAnalysis: не открыть " << filename << "\n";
        return;
    }

    out << "lambda_nm,R_cavity,T_cavity,R_periodic,T_periodic\n";
    out << std::scientific << std::setprecision(8);

    for (const auto& r : result_)
        out << r.lambda_nm  << ","
            << r.R_cavity   << ","
            << r.T_cavity   << ","
            << r.R_periodic << ","
            << r.T_periodic << "\n";

    std::cout << "CavityAnalysis: spectra written to " << filename << "\n";
}

void CavityAnalysis::writeFieldSnapshotsCSV(const std::string& prefix) const
{
    std::string infoFile = prefix + "_structure_info.csv";
    std::ofstream out(infoFile);
    if (!out.is_open()) return;

    out << "parameter,value\n";
    out << "structStart," << structStart_ << "\n";
    out << "structEnd,"   << structEnd_   << "\n";
    out << "cavStart,"    << cavStart_    << "\n";
    out << "cavEnd,"      << cavEnd_      << "\n";
    out << "lambda0_nm,"  << cfgC_.lambda0_nm << "\n";
    out << "nH,"          << cfgC_.nH     << "\n";
    out << "nL,"          << cfgC_.nL     << "\n";
    out << "nPairs,"      << cfgC_.nPairs << "\n";

    std::cout << "CavityAnalysis: structure info written to " << infoFile << "\n";
}


//#include <fstream>
//#include <iostream>
//#include <iomanip>
//#include <cmath>
//#include <algorithm>
//#include <stdexcept>
//#include <numeric>
//
//
//// Вычисляет толщины слоёв в ячейках
//// dH = round(lambda0 / (4*nH) / dx)
//// dL = round(lambda0 / (4*nL) / dx)
//// dCav = round(lambda0 / (2*nCav) / dx)
//void CavityAnalysis::computeLayerCells(int& dH, int& dL, int& dCav) const
//{
//    dH   = static_cast<int>(std::round(cfgC_.lambda0_nm / (4.0 * cfgC_.nH   * cfgC_.dx_nm)));
//    dL   = static_cast<int>(std::round(cfgC_.lambda0_nm / (4.0 * cfgC_.nL   * cfgC_.dx_nm)));
//
//    if (cfgC_.cavThickness_nm > 0.0)
//        dCav = static_cast<int>(std::round(cfgC_.cavThickness_nm / cfgC_.dx_nm));
//    else
//        dCav = static_cast<int>(std::round(cfgC_.lambda0_nm / (2.0 * cfgC_.nCav * cfgC_.dx_nm)));
//
//    std::cout << "CavityAnalysis geometry:\n"
//              << "  dH   = " << dH   << " cells (" << dH   * cfgC_.dx_nm << " nm)\n"
//              << "  dL   = " << dL   << " cells (" << dL   * cfgC_.dx_nm << " nm)\n"
//              << "  dCav = " << dCav << " cells (" << dCav * cfgC_.dx_nm << " nm)\n"
//              << "  nPairs = " << cfgC_.nPairs << "\n";
//}
//
//
//void CavityAnalysis::runVacuum(int monRefPos, int monTransPos,
//                                Monitor& mon,
//                                std::vector<double>& eTransIncEx,
//                                std::vector<double>& eTransIncHy)
//{
//    SimulationParameters pRef = base_;
//
//    FieldMonitor transMonVac(monTransPos);
//
//    FDTD1D simRef(pRef);
//    simRef.attachMonitor(&mon.incMonitor());
//    simRef.attachMonitor(&transMonVac);
//    simRef.run();
//
//    eTransIncEx = transMonVac.dataEx;
//    eTransIncHy = transMonVac.dataHy;
//
//    std::cout << "  Vacuum run done. incMon size=" << mon.incMonitor().size() << "\n";
//}
//
//// ДПФ: амплитуда |F(f)|
//static double dftAmp(const std::vector<double>& sig, double freq, double dt)
//{
//    double re = 0.0, im = 0.0;
//    const double omega = 2.0 * M_PI * freq;
//    const int N = static_cast<int>(sig.size());
//    for (int n = 0; n < N; ++n) {
//        const double ph = -omega * n * dt;
//        re += sig[n] * std::cos(ph);
//        im += sig[n] * std::sin(ph);
//    }
//    return std::sqrt(re * re + im * im);
//}
//
//
//// Один прогон с материалом
//void CavityAnalysis::runOnce(const MaterialLayout& layout,
//                              int monRefPos, int monTransPos,
//                              bool saveField,
//                              const std::vector<double>& eTransIncEx,
//                              const std::vector<double>& eTransIncHy,
//                              std::vector<double>& outR,
//                              std::vector<double>& outT,
//                              std::vector<double>& outFreqs)
//{
//    SimulationParameters pWork = base_;
//
//    Monitor      mon(pWork.dt, monRefPos);
//    FieldMonitor transMonWork(monTransPos);
//
//    {
//        SimulationParameters pRef = base_;
//        FDTD1D simRef(pRef);
//        simRef.attachMonitor(&mon.incMonitor());
//        simRef.run();
//    }
//
//    FDTD1D simWork(pWork, layout);
//    simWork.attachMonitor(&mon.totMonitor());
//    simWork.attachMonitor(&transMonWork);
//    simWork.run();
//
//
//    if (saveField) {
//        simWork.writeFieldCSV(cfgC_.fieldPrefix + "_raw.csv");
//        // Дополнительно пишем финальный срез Ex и Hy через мониторы
//    }
//
//    // Спектр R
//    const double fMin = cfgC_.dx_nm / cfgC_.lambdaMax;
//    const double fMax = cfgC_.dx_nm / cfgC_.lambdaMin;
//    const int    nF   = cfgC_.nWavelengths;
//
//    auto specR = mon.computeReflection(fMin, fMax, nF);
//
//    outFreqs.resize(nF);
//    outR.resize(nF);
//    outT.resize(nF);
//
//    const int Ntrans = static_cast<int>(
//        std::min({eTransIncEx.size(), eTransIncHy.size(),
//                  transMonWork.dataEx.size(), transMonWork.dataHy.size()}));
//
//
//
//    // T = |E_trans_work|^2 / |E_trans_inc|^2
//    std::vector<double> eTrW(transMonWork.dataEx.begin(),
//                                  transMonWork.dataEx.begin() + Ntrans);
//    std::vector<double> eTrI(eTransIncEx.begin(),
//                                  eTransIncEx.begin() + Ntrans);
//
//    for (int i = 0; i < nF; ++i) {
//        outFreqs[i] = specR[i].first;
//        outR[i]     = specR[i].second;
//    }
//
//    std::vector<double> incAmp(nF, 0.0);
//
//
//    double maxIncAmp = 0.0;
//    for (int i = 0; i < nF; ++i) {
//        const double f = outFreqs[i];
//        incAmp[i] = dftAmp(eTrI, f, pWork.dt);
//        if (incAmp[i] > maxIncAmp)
//            maxIncAmp = incAmp[i];
//    }
//
//    const double threshold = maxIncAmp * 1e-3;
//    // считаем T
//    for (int i = 0; i < nF; ++i) {
//        const double f    = outFreqs[i];
//        const double ampI = incAmp[i];
//        const double ampW = dftAmp(eTrW, f, pWork.dt);
//
//        if (ampI > threshold)
//            outT[i] = (ampW * ampW) / (ampI * ampI);
//        else
//            outT[i] = 0.0;
//    }
//}
//
//
//std::vector<CavityResult> CavityAnalysis::runCavity()
//{
//    int dH, dL, dCav;
//    computeLayerCells(dH, dL, dCav);
//
//    const int mirrorCells = cfgC_.nPairs * (dH + dL);
//    const int totalStruct = 2 * mirrorCells + dCav;
//
//    structStart_ = base_.source_pos + 200;
//    structEnd_   = structStart_ + totalStruct;
//    cavStart_    = structStart_ + mirrorCells;
//    cavEnd_      = cavStart_ + dCav;
//
//    if (structEnd_ >= base_.nx - base_.pmlThickness - 50)
//        throw std::runtime_error("CavityAnalysis: структура не помещается в сетку");
//
//    const int monRefPos   = base_.source_pos - cfgC_.monRefOffset;
//    const int monTransPos = structEnd_ + cfgC_.monTransOffset;
//
//    std::cout << "CavityAnalysis::run()\n"
//              << "  structStart=" << structStart_
//              << "  cavStart=" << cavStart_
//              << "  cavEnd=" << cavEnd_
//              << "  structEnd=" << structEnd_ << "\n"
//              << "  monRef=" << monRefPos
//              << "  monTrans=" << monTransPos << "\n";
//
//    std::vector<double> eTransIncEx, eTransIncHy;
//    {
//        Monitor monVac(base_.dt, monRefPos);
//        runVacuum(monRefPos, monTransPos, monVac, eTransIncEx, eTransIncHy);
//    }
//
//    // полная полость (DBR | cavity | DBR)
//    MaterialLayout layoutCav;
//    layoutCav.addBraggCavity(structStart_, cfgC_.nH, cfgC_.nL,
//                              dH, dL, cfgC_.nPairs, dCav);
//
//    std::vector<double> Rcav, Tcav, freqs;
//    runOnce(layoutCav, monRefPos, monTransPos, true,
//            eTransIncEx, eTransIncHy, Rcav, Tcav, freqs);
//
//    // чисто периодическая структура (без дефекта)
//    MaterialLayout layoutPer;
//    {
//        int cur = structStart_;
//        // 2*nPairs + 1 пара (без полости — просто (HL)^(2N+1))
//        const int totalPairs = 2 * cfgC_.nPairs + 1;
//        for (int p = 0; p < totalPairs; ++p) {
//            layoutPer.add(MaterialRegion(cur, cur + dH, cfgC_.nH * cfgC_.nH, 1.0, "H"));
//            cur += dH;
//            layoutPer.add(MaterialRegion(cur, cur + dL, cfgC_.nL * cfgC_.nL, 1.0, "L"));
//            cur += dL;
//        }
//    }
//
//    std::vector<double> Rper, Tper, freqsPer;
//    runOnce(layoutPer, monRefPos, monTransPos, false,
//            eTransIncEx, eTransIncHy, Rper, Tper, freqsPer);
//
//    const int nF = cfgC_.nWavelengths;
//    result_.clear();
//    result_.reserve(nF);
//
//    for (int i = nF - 1; i >= 0; --i) {
//        const double f   = freqs[i];
//        const double lam = (f > 1e-30) ? (cfgC_.dx_nm / f) : 0.0;
//
//        CavityResult r;
//        r.lambda_nm   = lam;
//        r.R_cavity    = Rcav[i];
//        r.T_cavity    = Tcav[i];
//        r.R_periodic  = (i < (int)Rper.size())  ? Rper[i]  : 0.0;
//        r.T_periodic  = (i < (int)Tper.size())  ? Tper[i]  : 0.0;
//        result_.push_back(r);
//    }
//
//    return result_;
//}
//
//void CavityAnalysis::writeSpectraCSV(const std::string& filename) const
//{
//    if (result_.empty()) {
//        std::cerr << "CavityAnalysis::writeSpectraCSV: нет данных\n";
//        return;
//    }
//    std::ofstream out(filename);
//    if (!out.is_open()) {
//        std::cerr << "CavityAnalysis: не открыть " << filename << "\n";
//        return;
//    }
//    out << "lambda_nm,R_cavity,T_cavity,R_periodic,T_periodic\n";
//    out << std::scientific << std::setprecision(8);
//    for (const auto& r : result_)
//        out << r.lambda_nm   << ","
//            << r.R_cavity    << ","
//            << r.T_cavity    << ","
//            << r.R_periodic  << ","
//            << r.T_periodic  << "\n";
//
//    std::cout << "CavityAnalysis: spectra written to " << filename << "\n";
//}
//
//void CavityAnalysis::writeFieldSnapshotsCSV(const std::string& prefix) const
//{
//    // Поле уже записано через simWork.writeFieldCSV() в runOnce()
//    // Здесь дополнительно пишем info-файл с координатами структуры
//    std::string infoFile = prefix + "_structure_info.csv";
//    std::ofstream out(infoFile);
//    if (!out.is_open()) return;
//
//    out << "parameter,value\n";
//    out << "structStart," << structStart_ << "\n";
//    out << "structEnd,"   << structEnd_   << "\n";
//    out << "cavStart,"    << cavStart_    << "\n";
//    out << "cavEnd,"      << cavEnd_      << "\n";
//    out << "lambda0_nm,"  << cfgC_.lambda0_nm << "\n";
//    out << "nH,"          << cfgC_.nH     << "\n";
//    out << "nL,"          << cfgC_.nL     << "\n";
//    out << "nPairs,"      << cfgC_.nPairs << "\n";
//
//    std::cout << "CavityAnalysis: structure info written to " << infoFile << "\n";
//}