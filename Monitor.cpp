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
    const int N = static_cast<int>(signal.size());
    std::vector<std::complex<double>> result(freqs.size());

    for (std::size_t fi = 0; fi < freqs.size(); ++fi) {
        const double omega = 2.0 * M_PI * freqs[fi];
        std::complex<double> acc(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            double phase = -omega * n * dt;
            acc += signal[n] * std::complex<double>(std::cos(phase), std::sin(phase));
        }
        result[fi] = acc;
    }
    return result;
}

std::vector<std::pair<double,double>> Monitor::computeReflection(double fMin,
                                                       double fMax, int nFreqs) const
{
    const auto& inc = incMon_.data;
    const auto& tot = totMon_.data;

    if (inc.empty() || tot.empty())
        throw std::runtime_error("Monitor: нет данных. Запустите оба прогона.");

    const int N = static_cast<int>(std::min(inc.size(), tot.size()));

    // E_ref[n] = E_tot[n] - E_inc[n]
    std::vector<double> eRef(N);
    for (int n = 0; n < N; ++n) {
        eRef[n] = tot[n] - inc[n];
        std::cout << eRef[n] << " = " << tot[n] << " - " << inc[n] << std::endl;
    }

    // Массив частот(равномерная секта)
    std::vector<double> freqs(nFreqs);
    const double df = (nFreqs > 1) ? (fMax - fMin) / (nFreqs - 1) : 0.0;
    for (int i = 0; i < nFreqs; ++i)
        freqs[i] = fMin + i * df;

    // ДПФ
    auto specInc = dft(inc,  dt_, freqs);
    auto specRef = dft(eRef, dt_, freqs);

    // R(f) = |E_ref(f)| / |E_inc(f)|
    double maxAmp = 0.0;
    for (auto& c : specInc) maxAmp = std::max(maxAmp, std::abs(c));
    const double threshold = maxAmp * 1e-3;

    std::vector<std::pair<double,double>> result(nFreqs);
    for (int i = 0; i < nFreqs; ++i) {
        double ampInc = std::abs(specInc[i]);
        double R = (ampInc > threshold) ? std::abs(specRef[i]) / ampInc : 0.0;
        result[i] = {freqs[i], R};
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

void Monitor::writeReflectionCSV(const std::string& filename,
                                  double fMin, double fMax, int nFreqs) const
{
    auto spec = computeReflection(fMin, fMax, nFreqs);

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Monitor: cannot open " << filename << "\n";
        return;
    }
    out << "freq,R\n";
    out << std::scientific;
    for (auto& [f, R] : spec)
        out << f << "," << R << "\n";

    std::cout << "Monitor: reflection spectrum written to " << filename << "\n";
}




void Monitor::writeTimeSeriesCSV(const std::string& filename) const
{
    const auto& inc = incMon_.data;
    const auto& tot = totMon_.data;

    if (inc.empty() || tot.empty()) {
        throw std::runtime_error("Monitor::writeTimeSeriesCSV: нет данных incident/total");
    }

    const int N = static_cast<int>(std::min(inc.size(), tot.size()));

    std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error("Monitor::writeTimeSeriesCSV: не удалось открыть файл " + filename);
    }

    out << "step,time,inc,tot,eref\n";
    out << std::scientific << std::setprecision(10);

    for (int n = 0; n < N; ++n) {
        const double t    = n * dt_;
        const double eRef = tot[n] - inc[n];
        out << n << "," << t << "," << inc[n] << "," << tot[n] << "," << eRef << "\n";
    }
}
