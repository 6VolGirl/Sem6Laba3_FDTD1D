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
