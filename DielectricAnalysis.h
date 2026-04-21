//
// Created by 6anna on 21.04.2026.
//

#ifndef DIELECTRICANALYSIS_H
#define DIELECTRICANALYSIS_H



#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include "SimulationParameters.h"
#include "Material.h"
#include "Monitor.h"
#include "FDTD1D.h"
#include "PMLAnalysis.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class DielectricAnalysis {
public:

    struct ConfigDiel {
        ConfigDiel() {}

        // Диапазон длин волн
        double lambdaMin  = 300.0;
        double lambdaMax  = 900.0;
        int    nWavelengths = 300;

        int monitorOffset = 100;


        double a_nm = 1.0;   // если dx = 1 нм, то a = 1 нм
    };
private:
    SimulationParameters base_;
    MaterialLayout      layout_;
    ConfigDiel          cfgD_;
    std::vector<std::pair<double,double>> result_;





public:

    DielectricAnalysis(const SimulationParameters& base, const MaterialLayout& layout,
                       ConfigDiel cfg = ConfigDiel())
        : base_(base), layout_(layout), cfgD_(cfg) {}


    // Возвращает вектор {lambda_nm, R}
    std::vector<std::pair<double,double>> run();

    // Теоретический R Френеля для нормального падения
    static double fresnelR(double n1, double n2) {
        double r = (n1 - n2) / (n1 + n2);
        return r * r;
    }

    // Теоретический R для всего диапазона
    double theoreticalR(double n1 = 1.0, double n2 = 1.45) const {
        return fresnelR(n1, n2);
    }


    void writeCSV(const std::string& filename,
                  double n1 = 1.0, double n2 = 1.45) const;


    const std::vector<std::pair<double,double>>& result() const {
        return result_;
    }

};




#endif //DIELECTRICANALYSIS_H
