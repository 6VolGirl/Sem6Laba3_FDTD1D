//
// Created by 6anna on 20.04.2026.
//

#ifndef FDTD1D_H
#define FDTD1D_H

#include <vector>
#include <string>
#include "SimulationParameters.h"
#include "Sources.h"
#include "PMLCoefficients.h"
#include "Monitor.h"

class FDTD1D {
    const SimulationParameters& p_;

    std::vector<double> Ex_;
    std::vector<double> Hy_;
    std::vector<double> eps_;
    std::vector<double> mu_;
    std::vector<double> sigmaE_;
    std::vector<double> sigmaM_;

    PMLSigma pml_;

    // Ex по времени для записи в CSV
    std::vector<std::vector<double>> snapshotsEx_;

    // Источники
    CWSource cw_;
    GaussianSource gauss_;

    double sourceValue(double t) const;

public:
    explicit FDTD1D(const SimulationParameters& p);

    void run();
    void writeFieldCSV(const std::string& filename) const;

};

#endif //FDTD1D_H


