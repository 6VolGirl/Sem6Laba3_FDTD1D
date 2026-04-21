
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
#include "PMLAnalysis.h"



static void runCase(const SimulationParameters& base,
                    SimulationParameters::SourceType st,
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
static void runCase2(const SimulationParameters& base,
                    SimulationParameters::SourceType st,
                    SimulationParameters::InjectionType it,
                    const std::string& filename,
                    const MaterialLayout& layout)
{
    SimulationParameters p = base;
    p.sourceType    = st;
    p.injectionType = it;

    FDTD1D sim(p, layout);
    sim.run();
    sim.writeFieldCSV(filename);
}

int main() {
    SimulationParameters base;

    // Сетка
    base.nx            = 1500;
    base.dx            = 1.0;
    base.courantNumber = 0.5;
    base.dt            = base.courantNumber * base.dx; // = 0.5
    base.numTimeSteps  = 2000;

    base.eps0 = 1.0;
    base.mu0  = 1.0;

    // Источник
    base.source_pos   = 200;
    //
    // const double fMin_src = 1.0 / 900.0;
    // const double fMax_src = 1.0 / 300.0;
    // base.sourceFreq   = 0.5 * (fMin_src + fMax_src);  // ~ 0.002222
    // base.sourceFWidth = fMax_src - fMin_src;
    //
    base.sourceFreq   = 0.05;   // период ~ 20 шагов по времени => ~40 ячеек на λ
    base.sourceFWidth = 0.02;
    base.sourceType     = SimulationParameters::GAUSS;
    base.injectionType  = SimulationParameters::SOFT;


    base.pmlThickness   = 20;
    base.pmlDamping     = 1e-4;
    base.pmlProfilePower = 3;


    base.snapshotEvery = 2;

    const int monitorPos = base.nx - base.pmlThickness - 15;

    try {
        //  Задания 0: запись полей для разных источников
        std::cout << "\n=== Tasks 1-2: field snapshots ===\n";
        runCase(base, SimulationParameters::CW,    SimulationParameters::SOFT,    "field_CW_soft.csv");
        runCase(base, SimulationParameters::CW,    SimulationParameters::CURRENT, "field_CW_current.csv");
        runCase(base, SimulationParameters::GAUSS, SimulationParameters::SOFT,    "field_Gauss_soft.csv");
        runCase(base, SimulationParameters::GAUSS, SimulationParameters::CURRENT, "field_Gauss_current.csv");


        //  Задания 1.3: анализ спектра отражения от PML
        std::cout << "\n=== Tasks 1.3: PML reflection analysis ===\n";


        PMLAnalysis::Config cfg;
        cfg.fMin          = 0.005;
        cfg.fMax          = 0.09;
        cfg.nFreqs        = 300;
        cfg.monitorOffset = 150;     // монитор на 150 ячеек левее источника

        PMLAnalysis analysis(base, cfg);

        // Задание 1.4: спектр R(f) для трёх профилей, width=20
        std::cout << "\n--- Spectra vs profile (width=20) ---\n";
        analysis.task4_spectraVsProfile(20, {0, 2, 3});

        // Дополнительно: спектр для нескольких толщин при кубическом профиле
        std::cout << "\n--- Spectra for cubic profile, varying width ---\n";
        for (int w : {5, 10, 20, 30}) {
            auto spec = analysis.runOne(w, 3);
            std::string fname = "reflection_cubic_width_" + std::to_string(w) + ".csv";

            std::ofstream out(fname);
            out << "freq,R\n" << std::scientific;
            for (auto& [f, R] : spec) out << f << "," << R << "\n";
            std::cout << fname << " written\n";
        }

        // Задание 1.5: R_max(width) для трёх профилей
        std::cout << "\n--- R_max vs width for profiles 0, 2, 3 ---\n";
        analysis.task5_maxRvsWidth({5, 8, 10, 15, 20, 25, 30}, {0, 2, 3});



//─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
        //  Задание 2: отражение от кварца R(λ)
        std::cout << "\n=== Task 2.3: Reflection from SiO2 half-space ===\n";


        // Материал: кварц в правой половине
        MaterialLayout layout;
        layout.add(MaterialRegion::Silica(base.nx / 2, base.nx));

        runCase2(base, SimulationParameters::GAUSS, SimulationParameters::SOFT, "field_Gauss_soft_Silica.csv", layout);



        std::cout << "All runs finished.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}