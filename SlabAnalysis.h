//
// Created by 6anna on 25.04.2026.
//

#ifndef SLABANALYSIS_H
#define SLABANALYSIS_H



#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include "SimulationParameters.h"
#include "Material.h"
#include "Monitor.h"
#include "FDTD1D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct SlabResult {
    double lambda_nm = 0.0;
    double R_fdtd    = 0.0;
    double T_fdtd    = 0.0;
    double R_theory  = 0.0;
    double T_theory  = 0.0;
};

class SlabAnalysis {
public:
    struct ConfigSlab {
        double lambdaMin    = 300.0;   // нм
        double lambdaMax    = 900.0;   // нм
        int    nWavelengths = 300;
        double a_nm         = 1.0;     // 1 ячейка = 1 нм
        int    monRefOffset = 100;     // монитор отражения отн источника
        int    monTransOffset = 100;   // монитор пропускания отн источника
        double n1 = 1.0;
        double n2 = 1.45;
    };

private:
    SimulationParameters base_;
    ConfigSlab cfgS_;
    std::vector<SlabResult>  result_;


public:
    explicit SlabAnalysis(const SimulationParameters& base)
            : base_(base), cfgS_() {}

    SlabAnalysis(const SimulationParameters& base, const ConfigSlab& cfgS)
        : base_(base), cfgS_(cfgS) {}


    // Запуск для одной толщины L (в нм)
    std::vector<SlabResult> runSlab(double L_nm, int numSteps = 0);

    // Запуск для нескольких L
    void runMultipleL(const std::vector<double>& L_vals,
                      const std::string& prefix = "slab_");

    // Теоретические R и T
    static double theoreticalR(double lambda_nm, double L_nm,
                                double n1, double n2);
    static double theoreticalT(double lambda_nm, double L_nm,
                                double n1, double n2);

    static void writeCSV(const std::string& filename,
                         const std::vector<SlabResult>& data);

    const std::vector<SlabResult>& result() const {
        return result_;
    }

};



#endif //SLABANALYSIS_H
