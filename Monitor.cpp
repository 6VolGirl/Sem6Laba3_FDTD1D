//
// Created by 6anna on 20.04.2026.
//

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include "Monitor.h"


//ДПФ: E(f) = sum_n  signal[n] * exp(-i 2π f n dt)

std::vector<std::complex<double>> Monitor::dft(const std::vector<double>& signal, double dt,
             const std::vector<double>& freqs)
{
    const int N = (int)signal.size();
    std::vector<std::complex<double>> result(freqs.size());

    for (std::size_t fi = 0; fi < freqs.size(); ++fi) {
        const double omega = 2.0 * M_PI * freqs[fi];
        std::complex<double> acc(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            double phase = -omega * n * dt;
            acc += signal[n] * std::complex<double>(std::cos(phase),
                                                    std::sin(phase));
        }
        result[fi] = acc;
    }
    return result;
}


std::vector<std::pair<double,double>>
Monitor::computeReflection(double fMin, double fMax, int nFreqs) const
{
    const int N = std::min(incMon_.size(), totMon_.size());
    if (N == 0)
        throw std::runtime_error("Monitor::computeReflection: нет данных.");


    // S нелинеен — вычитание S_tot - S_inc некорректно для смешанных волн
    std::vector<double> eInc(N), eRef(N);
    for (int n = 0; n < N; ++n) {
        eInc[n] = incMon_.dataEx[n];
        eRef[n] = totMon_.dataEx[n] - incMon_.dataEx[n];  // E_ref = E_tot - E_inc
    }

    std::vector<double> freqs(nFreqs);
    const double df = (nFreqs > 1) ? (fMax - fMin) / (nFreqs - 1) : 0.0;
    for (int i = 0; i < nFreqs; ++i)
        freqs[i] = fMin + i * df;

    auto specInc = dft(eInc, dt_, freqs);
    auto specRef = dft(eRef, dt_, freqs);

    double maxAmpInc = 0.0;
    for (auto& c : specInc)
        maxAmpInc = std::max(maxAmpInc, std::abs(c));
    const double threshold = maxAmpInc * 1e-3;

    std::vector<std::pair<double,double>> result(nFreqs);
    for (int i = 0; i < nFreqs; ++i) {
        double ampInc = std::abs(specInc[i]);
        double ampRef = std::abs(specRef[i]);
        // R = |E_ref|^2 / |E_inc|^2
        double R = (ampInc > threshold)
                   ? (ampRef * ampRef) / (ampInc * ampInc)
                   : 0.0;
        result[i] = {freqs[i], R};
    }
    return result;
}

std::vector<std::pair<double,double>>
Monitor::computeTransmission(const FieldMonitor& transWork,
                             const FieldMonitor& transInc,
                             double fMin, double fMax, int nFreqs) const
{
    const int N = std::min({incMon_.size(),
                            transWork.size(),
                            transInc.size()});
    if (N == 0)
        throw std::runtime_error("Monitor::computeTransmission: нет данных.");

    std::vector<double> sInc(N), sTrans(N);
    for (int n = 0; n < N; ++n) {
        sInc[n]   = transInc.dataEx[n]  * transInc.dataHy[n];
        sTrans[n] = transWork.dataEx[n] * transWork.dataHy[n];
    }

    std::vector<double> freqs(nFreqs);
    const double df = (nFreqs > 1) ? (fMax - fMin) / (nFreqs - 1) : 0.0;
    for (int i = 0; i < nFreqs; ++i)
        freqs[i] = fMin + i * df;

    auto specInc   = dft(sInc,   dt_, freqs);
    auto specTrans = dft(sTrans, dt_, freqs);

    double maxAmpInc = 0.0;
    for (auto& c : specInc)
        maxAmpInc = std::max(maxAmpInc, std::abs(c));
    const double threshold = maxAmpInc * 1e-3;

    std::vector<std::pair<double,double>> result(nFreqs);
    for (int i = 0; i < nFreqs; ++i) {
        double ampInc = std::abs(specInc[i]);
        double T = (ampInc > threshold)
                   ? std::abs(specTrans[i]) / ampInc
                   : 0.0;
        result[i] = {freqs[i], T};
    }
    return result;
}



std::vector<std::pair<double,double>>
Monitor::computeReflectionPoynting(double fMin, double fMax, int nFreqs) const
{
    const int N = std::min(incMon_.size(), totMon_.size());
    if (N == 0)
        throw std::runtime_error("Monitor::computeReflectionPoynting: нет данных.");

    // Разностные поля: E_ref = E_tot - E_inc
    std::vector<double> exRef(N), hyRef(N);
    for (int n = 0; n < N; ++n) {
        exRef[n] = totMon_.dataEx[n] - incMon_.dataEx[n];
        hyRef[n] = totMon_.dataHy[n] - incMon_.dataHy[n];
    }

    std::vector<double> freqs(nFreqs);
    const double df = (nFreqs > 1) ? (fMax - fMin) / (nFreqs - 1) : 0.0;
    for (int i = 0; i < nFreqs; ++i)
        freqs[i] = fMin + i * df;

    const std::vector<double> exIncVec(incMon_.dataEx.begin(), incMon_.dataEx.begin() + N);
    const std::vector<double> hyIncVec(incMon_.dataHy.begin(), incMon_.dataHy.begin() + N);

    auto specExInc = dft(exIncVec, dt_, freqs);
    auto specHyInc = dft(hyIncVec, dt_, freqs);
    auto specExRef = dft(exRef,    dt_, freqs);
    auto specHyRef = dft(hyRef,    dt_, freqs);

    // ── ИСПРАВЛЕНИЕ: порог через |S_inc|, а не просто S_inc > 0
    double maxSinc = 0.0;
    for (int i = 0; i < nFreqs; ++i) {
        double s = std::abs(0.5 * std::real(specExInc[i] * std::conj(specHyInc[i])));
        if (s > maxSinc) maxSinc = s;
    }
    const double threshold = maxSinc * 1e-6;

    std::vector<std::pair<double,double>> result(nFreqs);
    for (int i = 0; i < nFreqs; ++i) {
        const double S_inc = 0.5 * std::real(specExInc[i] * std::conj(specHyInc[i]));
        const double S_ref = 0.5 * std::real(specExRef[i] * std::conj(specHyRef[i]));

        double R = 0.0;
        if (std::abs(S_inc) > threshold) {
            // ── ИСПРАВЛЕНИЕ: R = -S_ref / S_inc
            // S_ref < 0 (волна идёт влево), S_inc > 0 (волна идёт вправо)
            // поэтому знак "-" даёт R >= 0
            R = -S_ref / S_inc;
            if (R < 0.0) R = 0.0; // подавить численный шум
        }
        result[i] = {freqs[i], R};
    }
    return result;
}

std::vector<std::pair<double,double>>
Monitor::computeTransmissionPoynting(const FieldMonitor& transWork,
                                     const FieldMonitor& transInc,
                                     double fMin, double fMax, int nFreqs) const
{
    const int N = std::min(transWork.size(), transInc.size());
    if (N == 0)
        throw std::runtime_error("Monitor::computeTransmissionPoynting: нет данных.");

    const std::vector<double> exWork(transWork.dataEx.begin(),
                                     transWork.dataEx.begin() + N);
    const std::vector<double> hyWork(transWork.dataHy.begin(),
                                     transWork.dataHy.begin() + N);
    const std::vector<double> exInc(transInc.dataEx.begin(),
                                    transInc.dataEx.begin() + N);
    const std::vector<double> hyInc(transInc.dataHy.begin(),
                                    transInc.dataHy.begin() + N);

    // Сетка частот
    std::vector<double> freqs(nFreqs);
    const double df = (nFreqs > 1) ? (fMax - fMin) / (nFreqs - 1) : 0.0;
    for (int i = 0; i < nFreqs; ++i)
        freqs[i] = fMin + i * df;

    // ДПФ: Ex и Hy отдельно для обоих мониторов
    auto specExWork = dft(exWork, dt_, freqs);
    auto specHyWork = dft(hyWork, dt_, freqs);
    auto specExInc  = dft(exInc,  dt_, freqs);
    auto specHyInc  = dft(hyInc,  dt_, freqs);

    // Порог по максимуму S_inc_trans
    double maxSinc = 0.0;
    for (int i = 0; i < nFreqs; ++i) {
        double s = 0.5 * std::real(specExInc[i] * std::conj(specHyInc[i]));
        if (s > maxSinc) maxSinc = s;
    }
    const double threshold = maxSinc * 1e-3;

    std::vector<std::pair<double,double>> result(nFreqs);
    for (int i = 0; i < nFreqs; ++i) {
        const double S_inc_trans = 0.5 * std::real(specExInc[i]  * std::conj(specHyInc[i]));
        const double S_trans     = 0.5 * std::real(specExWork[i] * std::conj(specHyWork[i]));

        double T = 0.0;
        if (S_inc_trans > threshold)
            T = S_trans / S_inc_trans;
        if (T < 0.0) T = 0.0;   // численный шум

        result[i] = {freqs[i], T};
    }
    return result;
}


double Monitor::maxReflection(double fMin, double fMax, int nFreqs) const
{
    auto spec = computeReflection(fMin, fMax, nFreqs);
    double rMax = 0.0;
    for (auto& [f, R] : spec)
        rMax = std::max(rMax, R);
    return rMax;
}


void Monitor::writeTimeSeriesCSV(const std::string& filename) const
{
    const int N = std::min(incMon_.size(), totMon_.size());
    if (N == 0)
        throw std::runtime_error("Monitor::writeTimeSeriesCSV: нет данных.");

    std::ofstream out(filename);
    if (!out.is_open())
        throw std::runtime_error("Monitor::writeTimeSeriesCSV: не открыть " + filename);

    out << "step,time,"
        << "Ex_inc,Hy_inc,S_inc,"
        << "Ex_tot,Hy_tot,S_tot,"
        << "Ex_ref,Hy_ref,S_ref\n";
    out << std::scientific << std::setprecision(10);

    for (int n = 0; n < N; ++n) {
        const double t      = n * dt_;
        const double exInc  = incMon_.dataEx[n];
        const double hyInc  = incMon_.dataHy[n];
        const double sInc   = exInc * hyInc;
        const double exTot  = totMon_.dataEx[n];
        const double hyTot  = totMon_.dataHy[n];
        const double sTot   = exTot * hyTot;
        const double exRef  = exTot - exInc;
        const double hyRef  = hyTot - hyInc;
        const double sRef   = exRef * hyRef;

        out << n      << "," << t      << ","
            << exInc  << "," << hyInc  << "," << sInc  << ","
            << exTot  << "," << hyTot  << "," << sTot  << ","
            << exRef  << "," << hyRef  << "," << sRef  << "\n";
    }
    std::cout << "Monitor: timeseries written to " << filename << "\n";
}


void Monitor::writeSlabTimeSeriesCSV(const std::string& filename,
                                      const FieldMonitor& transMon) const
{
    const int N = std::min({incMon_.size(), totMon_.size(), transMon.size()});
    if (N == 0) return;

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Monitor::writeSlabTimeSeriesCSV: не открыть " << filename << "\n";
        return;
    }

    out << "step,time,"
        << "Ex_inc,Hy_inc,S_inc,"
        << "Ex_tot,Hy_tot,S_tot,"
        << "Ex_ref,Hy_ref,S_ref,"
        << "Ex_trans,Hy_trans,S_trans\n";
    out << std::scientific << std::setprecision(10);

    for (int n = 0; n < N; ++n) {
        const double t      = n * dt_;
        const double exInc  = incMon_.dataEx[n];
        const double hyInc  = incMon_.dataHy[n];
        const double sInc   = exInc * hyInc;
        const double exTot  = totMon_.dataEx[n];
        const double hyTot  = totMon_.dataHy[n];
        const double sTot   = exTot * hyTot;
        const double exRef  = exTot - exInc;
        const double hyRef  = hyTot - hyInc;
        const double sRef   = exRef * hyRef;
        const double exTr   = transMon.dataEx[n];
        const double hyTr   = transMon.dataHy[n];
        const double sTr    = exTr * hyTr;

        out << n     << "," << t     << ","
            << exInc << "," << hyInc << "," << sInc << ","
            << exTot << "," << hyTot << "," << sTot << ","
            << exRef << "," << hyRef << "," << sRef << ","
            << exTr  << "," << hyTr  << "," << sTr  << "\n";
    }
    std::cout << "Monitor: slab timeseries written to " << filename << "\n";
}
