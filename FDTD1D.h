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

    // Поля на Yee-сетке
    std::vector<double> Ex_;   // размер nx+1
    std::vector<double> Hy_;   // размер nx   (полушаг)

    // Снимки Ex по времени для записи в CSV
    std::vector<std::vector<double>> snapshotsEx_;

    // Источники
    CWSource       cw_;
    GaussianSource gauss_;

    double sourceValue(double t) const;

public:
    explicit FDTD1D(const SimulationParameters& p);

    // Запуск главного временного цикла
    void run();

    // Запись массива Ex(x, t) в CSV c заголовком time_over_fL,x_tilde,Ez
    void writeFieldCSV(const std::string& filename) const;
};

#endif //FDTD1D_H


