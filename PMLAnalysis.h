//
// Created by 6anna on 21.04.2026.
//

#ifndef PMLANALYSIS_H
#define PMLANALYSIS_H


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

#include "SimulationParameters.h"
#include "FDTD1D.h"
#include "Monitor.h"


class PMLAnalysis {
public:
    // Параметры анализа
    struct Config {
        Config() {}
        double fMin   = 0.01;
        double fMax   = 0.09;
        int    nFreqs = 200;    // число точек по частоте
        int    monitorOffset = 50; // смещение монитора от источника влево
    };


private:
    SimulationParameters base_;
    Config cfg_;

    static void writePairCSV(const std::string& fname,
                             const std::string& header,
                             const std::vector<std::pair<double,double>>& data);


public:

    explicit PMLAnalysis(const SimulationParameters& base, Config cfg = {})
        : base_(base), cfg_(cfg) {}

    // Одна пара прогонов для заданного профиля и толщины
    std::vector<std::pair<double,double>> runOne(int pmlThickness, int profilePower) const;

    // спектр R(f) для трёх профилей и fix толщины
   void task4_spectraVsProfile(int pmlWidth,
                                const std::vector<int>& profiles) const;

    //  R_max(width) для нескольких профилей
    void task5_maxRvsWidth(const std::vector<int>& widths,
                           const std::vector<int>& profiles) const;


};



#endif //PMLANALYSIS_H
