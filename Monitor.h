//
// Created by 6anna on 20.04.2026.
//

#ifndef MONITOR_H
#define MONITOR_H

#include <vector>
#include <complex>
#include <string>
#include <cmath>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


class FieldMonitor {
public:
    int pos;
    std::vector<double> dataEx;  // Ex[pos][n]
    std::vector<double> dataHy;

    explicit FieldMonitor(int position) : pos(position) {}

    void record(const std::vector<double>& Ex, const std::vector<double>& Hy) {
        const int i = pos;
        const double hy = (i > 0 && i < (int)Hy.size())
                          ? 0.5 * (Hy[i - 1] + Hy[i])
                          : Hy[std::min(i, (int)Hy.size() - 1)];
        dataEx.push_back(Ex[i]);
        dataHy.push_back(hy);
    }

    std::vector<double> poynting() const {
        const int N = (int)std::min(dataEx.size(), dataHy.size());
        std::vector<double> s(N);
        for (int n = 0; n < N; ++n)
            s[n] = dataEx[n] * dataHy[n];
        return s;
    }

    void clear() { dataEx.clear(); dataHy.clear(); }

    int size() const { return (int)dataEx.size(); }
};



class Monitor {
private:
    double dt_;
    int    monPos_;

    FieldMonitor incMon_;
    FieldMonitor totMon_;


    // Возвращает спектр для freqs
    static std::vector<std::complex<double>>
    dft(const std::vector<double>& signal, double dt,
        const std::vector<double>& freqs);


public:
    Monitor(double dt, int monPos)
        : dt_(dt), monPos_(monPos), incMon_(monPos), totMon_(monPos) {}


    // incMon — прогон без PML
    // totMon — прогон
    FieldMonitor& incMonitor()  { return incMon_; }
    FieldMonitor& totMonitor()  { return totMon_; }

    int monitorPos() const { return monPos_; }

    // Вычисляет спектр R(f) после двух прогонов
    std::vector<std::pair<double,double>> computeReflection(double fMin,
                                           double fMax, int nFreqs) const;

    // Вычисляет спектр T(f) после двух прогонов
    std::vector<std::pair<double,double>> computeTransmission(
    const FieldMonitor& incMon,
    const FieldMonitor& transMon,
    double fMin, double fMax, int nFreqs) const;

    // Максимальное R в заданном диапазоне частот
    double maxReflection(double fMin, double fMax, int nFreqs) const;


    // Временные ряды: Ex, Hy, S = Ex*Hy
    void writeTimeSeriesCSV(const std::string& filename) const;

    // Запись для пластины с trans-монитором
    void writeSlabTimeSeriesCSV(const std::string& filename,
                                const FieldMonitor& transMon) const;


};

#endif // MONITOR_H



