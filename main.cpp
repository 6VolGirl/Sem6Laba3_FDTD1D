
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
#include "DielectricAnalysis.h"
#include "SlabAnalysis.h"



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

int main() {
    SimulationParameters base;

    // Сетка
    base.nx            = 3000;
    base.dx            = 1.0;
    base.courantNumber = 0.5;
    base.dt            = base.courantNumber * base.dx;
    base.numTimeSteps  = 5000;

    base.eps0 = 1.0;
    base.mu0  = 1.0;

    // Источник
    base.source_pos   = 1000;
    base.sourceFreq   = 0.05;
    base.sourceFWidth = 0.02;
    base.sourceType     = SimulationParameters::GAUSS;
    base.injectionType  = SimulationParameters::CURRENT;


    base.pmlThickness   = 50;
    base.pmlDamping     = 1e-12;
    base.pmlProfilePower = 3;


    base.snapshotEvery = 2;


    try {

        // //  Задания 0: запись полей для разных источников
        // {        std::cout << "\n=== Tasks 1-2: field snapshots ===\n";
        //     //runCase(base, SimulationParameters::CW,    SimulationParameters::SOFT,    "field_CW_soft.csv");
        //     //runCase(base, SimulationParameters::CW,    SimulationParameters::CURRENT, "field_CW_current.csv");
        //     runCase(base, SimulationParameters::GAUSS, SimulationParameters::SOFT,    "field_Gauss_soft.csv");
        //     //runCase(base, SimulationParameters::GAUSS, SimulationParameters::CURRENT, "field_Gauss_current.csv");
        //
        //     // Материал: кварц в правой половине
        //     MaterialLayout layout;
        //     layout.add(MaterialRegion::Silica(base.nx / 2, base.nx));
        //
        //     SimulationParameters p = base;
        //     FDTD1D sim(p, layout);
        //     sim.run();
        //     sim.writeFieldCSV("field_Gauss_soft_PML.csv");
        // }

        // //  Задания 1: анализ спектра отражения от PML
        // {std::cout << "\n=== Tasks 3: PML reflection analysis ===\n";
        //
        //     PMLAnalysis::Config cfg;
        //     cfg.fMin          = 0.005;
        //     cfg.fMax          = 0.4;
        //     cfg.nFreqs        = 300;
        //     cfg.monitorOffset = 150;     // монитор на 150 ячеек левее источника
        //
        //     PMLAnalysis analysis(base, cfg);
        //
        //     // Задание 1.4: спектр R(f) для трёх профилей, width=20
        //     std::cout << "\n--- Spectra vs profile (width=20) ---\n";
        //     analysis.task4_spectraVsProfile(20, {0, 2, 3});
        //
        //     std::cout << "\n--- Spectra for cubic profile, varying width ---\n";
        //     for (int w : {5, 10, 20, 30}) {
        //         auto spec = analysis.runOne(w, 3);
        //         std::string fname = "reflection_cubic_width_" + std::to_string(w) + ".csv";
        //
        //         std::ofstream out(fname);
        //         out << "freq,R\n" << std::scientific;
        //         for (auto& [f, R] : spec) out << f << "," << R << "\n";
        //         std::cout << fname << " written\n";
        //     }
        //
        //     // Задание 1.5: R_max(width) для трёх профилей
        //     std::cout << "\n--- R_max vs width for profiles 0, 2, 3 ---\n";
        //     analysis.task5_maxRvsWidth({5, 8, 10, 15, 20, 25, 30}, {0, 2, 3});
        // }

        // //  Задания 2: анализ спектра отражения от кварца (n = 1.45)
        // {
        //
        //     SimulationParameters pDiel = base;
        //     // Материал: кварц в правой половине
        //     MaterialLayout layout;
        //     layout.add(MaterialRegion::Silica(pDiel.nx / 2, pDiel.nx));
        //
        //     DielectricAnalysis::ConfigDiel cfgD;
        //     cfgD.lambdaMin     = 300.0;    // нм
        //     cfgD.lambdaMax     = 900.0;    // нм
        //     cfgD.nWavelengths  = 300;
        //     cfgD.monitorOffset = 100;
        //     cfgD.a_nm          = 1.0;
        //
        //     pDiel.nx = 2000;
        //     pDiel.source_pos = 600;
        //     pDiel.numTimeSteps = 10000;
        //     pDiel.sourceFreq   = 0.5 * (1.0/cfgD.lambdaMax + 1.0/cfgD.lambdaMin);
        //     pDiel.sourceFWidth = 1.0/cfgD.lambdaMin - 1.0/cfgD.lambdaMax;
        //
        //     DielectricAnalysis analysisDiel(base, layout, cfgD);
        //     analysisDiel.run();
        //
        //     // SimulationParameters p = base;
        //     // FDTD1D sim(p, layout);
        //     // sim.run();
        //     // sim.writeFieldCSV("field_Gauss_soft_PML.csv");
        //     analysisDiel.writeCSV("silica_reflection.csv", 1.0,1.45);
        //
        // }

        //  Задания 3: анализ спектра отражения и пропускания кварцевой пластины (n = 1.45)
        {
            SimulationParameters pSlab = base;



            SlabAnalysis::ConfigSlab cfgSlab;
            cfgSlab.lambdaMin     = 300.0;
            cfgSlab.lambdaMax     = 900.0;
            cfgSlab.nWavelengths  = 300;
            cfgSlab.a_nm          = 1.0;
            cfgSlab.monRefOffset  = 100;  // монитор R: source_pos - 100
            cfgSlab.monTransOffset = 100; // монитор T: slabEnd + 50
            cfgSlab.n1 = 1.0;
            cfgSlab.n2 = 1.45;

            pSlab.sourceFreq   = 0.5 * (1.0/cfgSlab.lambdaMax + 1.0/cfgSlab.lambdaMin);
            pSlab.sourceFWidth = 1.0/cfgSlab.lambdaMin - 1.0/cfgSlab.lambdaMax;
            pSlab.nx            = 2000;
            pSlab.source_pos    = 400;
            pSlab.numTimeSteps  = 16000;

            SlabAnalysis slabAn(pSlab, cfgSlab);

            // Задание 3.4: разные L
            slabAn.runMultipleL({100.0, 200.0, 400.0, 800.0}, "slab_");
        }

        std::cout << "All runs finished.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}