
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <string>
#include "SimulationParameters.h"
#include "FDTD1D.h"

static void runCase(const SimulationParameters& base,
                    SimulationParameters::SourceType   st,
                    SimulationParameters::InjectionType it,
                    const std::string& filename)
{
    SimulationParameters p = base;
    p.sourceType    = st;
    p.injectionType = it;

    FDTD1D sim(p);
    sim.run();
    sim.writeFieldCSV(filename);
}

int main() {
    SimulationParameters base;

    // Сетка
    base.nx            = 1000;
    base.dx            = 1.0;
    base.courantNumber = 0.5;
    base.dt            = base.courantNumber * base.dx; // = 0.5
    base.numTimeSteps  = 2000;

    base.eps0 = 1.0;
    base.mu0  = 1.0;

    // Источник
    base.source_pos   = base.nx / 2;
    base.sourceFreq   = 0.05;   // период ~ 20 шагов по времени => ~40 ячеек на λ
    base.sourceFWidth = 0.02;


    base.snapshotEvery = 2;

    try {
        // 1) CW, мягкий источник
        runCase(base, SimulationParameters::CW,    SimulationParameters::SOFT,
                "field_CW_soft.csv");

        // 2) CW, через плотность тока J
        runCase(base, SimulationParameters::CW,    SimulationParameters::CURRENT,
                "field_CW_current.csv");

        // 3) Гаусс, мягкий источник
        runCase(base, SimulationParameters::GAUSS, SimulationParameters::SOFT,
                "field_Gauss_soft.csv");

        // 4) Гаусс, через плотность тока J
        runCase(base, SimulationParameters::GAUSS, SimulationParameters::CURRENT,
                "field_Gauss_current.csv");

        std::cout << "All 4 runs finished. CSV format: time_over_fL,x_tilde,Ez\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}