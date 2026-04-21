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
                             const std::vector<std::pair<double,double>>& data)
    {
        std::ofstream out(fname);
        if (!out.is_open()) { std::cerr << "Cannot open " << fname << "\n"; return; }
        out << header << "\n" << std::scientific;
        for (auto& [a, b] : data) out << a << "," << b << "\n";
        std::cout << fname << " written\n";
    }



public:

    explicit PMLAnalysis(const SimulationParameters& base, Config cfg = {})
        : base_(base), cfg_(cfg) {}

    // Одна пара прогонов для заданного профиля и толщины
    std::vector<std::pair<double,double>>
    runOne(int pmlThickness, int profilePower) const
    {
        //  (задержанное на время 2*(dist_src_to_rightPML)/c)????!!!!!
        const int monPos = base_.source_pos - cfg_.monitorOffset;

        Monitor mon(base_.dt, monPos);

        SimulationParameters pRef = base_;
        pRef.pmlThickness    = 30;      // достаточно широкий
        pRef.pmlDamping      = 1e-14;   // "идеальное" поглощение
        pRef.pmlProfilePower = 3;

        FDTD1D sim(pRef);
        sim.attachMonitor(&mon.incMonitor());
        sim.run();




        // E_tot(t) = E_inc(t) + E_ref(t)
        {
            SimulationParameters pWork = base_;
            pWork.pmlThickness    = pmlThickness;
            pWork.pmlProfilePower = profilePower;
            // pmlDamping берётся из base_

            FDTD1D sim(pWork);
            sim.attachMonitor(&mon.totMonitor());
            sim.run();
        }

        return mon.computeReflection(cfg_.fMin, cfg_.fMax, cfg_.nFreqs);
    }

    // спектр R(f) для трёх профилей и fix толщины
   void task4_spectraVsProfile(int pmlWidth,
                                const std::vector<int>& profiles) const
    {
        for (int m : profiles) {
            auto spec = runOne(pmlWidth, m);

            std::string fname = "reflection_profile_" + std::to_string(m)
                              + "_width_" + std::to_string(pmlWidth) + ".csv";
            writePairCSV(fname, "freq,R", spec);
        }
    }

    //  R_max(width) для нескольких профилей
    void task5_maxRvsWidth(const std::vector<int>& widths,
                           const std::vector<int>& profiles) const
    {
        std::ofstream out("reflection_maxR_vs_width.csv");
        if (!out.is_open()) {
            std::cerr << "Cannot open reflection_maxR_vs_width.csv\n";
            return;
        }

        // Заголовок
        out << "width";
        for (int m : profiles) out << ",R_profile" << m;
        out << "\n";

        out << std::scientific;

        for (int w : widths) {
            out << w;
            for (int m : profiles) {
                auto spec = runOne(w, m);
                double rMax = 0.0;
                for (auto& [f, R] : spec) rMax = std::max(rMax, R);
                out << "," << rMax;
                std::cout << "  width=" << w << " profile=" << m
                          << " Rmax=" << rMax << "\n";
            }
            out << "\n";
        }
        std::cout << "reflection_maxR_vs_width.csv written\n";
    }


};



#endif //PMLANALYSIS_H
