//
// Created by 6anna on 26.04.2026.
//

#ifndef PHOTONCRYSTALANALYSIS_H
#define PHOTONCRYSTALANALYSIS_H


#include <vector>
#include <string>
#include <cmath>
#include "SimulationParameters.h"
#include "Material.h"
#include "Monitor.h"
#include "FDTD1D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct PCResult {
    double lambda_nm;
    double R_fdtd;
    double T_fdtd;
    double balance;
};

class PhotonicCrystalAnalysis {
public:

    struct ConfigPhot {
        // Оптические параметры
        double n_A = 2.28;     // TiO2
        double n_B = 1.45;     // SiO2

        // Толщины слоёв в нм
        double d_A_equal  = 100.0;   // TiO2, равные физ. толщины
        double d_B_equal  = 100.0;   // SiO2, равные физ. толщины
        double d_A_opt    = 55.0;    // TiO2, равные оптические пути
        double d_B_opt    = 100.0;   // SiO2, равные оптические пути

        double lambdaMin  = 300.0;
        double lambdaMax  = 900.0;
        int    nWavelengths = 400;

        double a_nm = 1.0;           // 1 ячейка = 1 нм

        int monRefOffset   = 150;    // левый монитор (R)
        int monTransOffset = 50;     // правый монитор (T)

        // вакуум с обеих сторон
        double n_outside = 1.0;
    };

private:
    const SimulationParameters& base_;
    ConfigPhot cfg_;

    // Вспомогательный ДПФ
    static double dftAmp(const std::vector<double>& sig, int N,
                         double freq, double dt);

    // Построить MaterialLayout из N периодов (A/B)
    static MaterialLayout buildLayout(int crystalStart,
                                       int dA_cells, int dB_cells,
                                       int numPeriods);


    std::vector<PCResult> runOnce(int crystalStart,
                                   int dA_cells, int dB_cells,
                                   int numPeriods) const;

public:
    PhotonicCrystalAnalysis(const SimulationParameters& base);
    PhotonicCrystalAnalysis(const SimulationParameters& base, const ConfigPhot& cfg);

    // Задание 4.2: равные физические толщины
    void runEqualThickness(const std::vector<int>& periodCounts,
                            const std::string& outPrefix);

    // Задание 4.2: равные оптические пути
    void runEqualOptPath(const std::vector<int>& periodCounts,
                          const std::string& outPrefix);

    // Записать один результат в CSV
    static void writeCSV(const std::string& filename,
                          const std::vector<PCResult>& data);
};


#endif //PHOTONCRYSTALANALYSIS_H
