//
// Created by 6anna on 20.04.2026.
//

#ifndef MONITOR_H
#define MONITOR_H

#include <vector>
#include <complex>
#include <string>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


class FieldMonitor {
public:
    int pos;
    std::vector<double> data;  // Ex[pos][n]

    explicit FieldMonitor(int position) : pos(position) {}

    void record(const std::vector<double>& Ex) {
        data.push_back(Ex[pos]);
    }

    void clear() { data.clear(); }
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

    // Максимальное R в заданном диапазоне частот
    double maxReflection(double fMin, double fMax, int nFreqs) const;


    void writeReflectionCSV(const std::string& filename,
                            double fMin, double fMax, int nFreqs) const;

    // Временные ряды
    void writeTimeSeriesCSV(const std::string& filename) const;


};

#endif // MONITOR_H



