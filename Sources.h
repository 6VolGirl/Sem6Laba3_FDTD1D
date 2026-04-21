//
// Created by 6anna on 20.04.2026.
//

#ifndef SOURCES_H
#define SOURCES_H



#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Непрерывная волна: sin(2π f t)
struct CWSource {
    double freq;
    explicit CWSource(double f) : freq(f) {}
    double operator()(double t) const {
        return std::sin(2.0 * M_PI * freq * t);
    }
};

// Гауссов импульс
struct GaussianSource {
    double freq;
    double fwidth;     // ширина спектра
    double w;          // длительность огибающей
    double t0;         // центр импульса
    double start_time;
    double finish_time;
    double cutoff;

    GaussianSource(double freq_, double fwidth_, double start = 0.0, double cutoff_ = 4.0)
    : freq(freq_), fwidth(fwidth_), start_time(start), cutoff(cutoff_) {
        w = 1.0 / fwidth;
        t0 = start_time + cutoff * w;
        finish_time = start_time + 2.0 * cutoff * w;
    }

    double operator()(double t) const {
        if (t < start_time || t > finish_time) return 0.0;
        double tau = t - t0;
        double env = std::exp(-0.5 * (tau * tau) / (w * w));
        return env * std::sin(2.0 * M_PI * freq * tau);
    }
};



#endif //SOURCES_H
