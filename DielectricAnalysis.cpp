//
// Created by 6anna on 24.04.2026.
//

#include "DielectricAnalysis.h"


#include <fstream>
#include <iostream>
#include <stdexcept>



std::vector<std::pair<double,double>> DielectricAnalysis::run()
{
    const int monPos = base_.source_pos - cfgD_.monitorOffset;
    if (monPos <= 0)
        throw std::runtime_error(
            "DielectricAnalysis: monitorOffset slishkom bolshoy");

    Monitor mon(base_.dt, monPos);

    // нет диэлектрика
    SimulationParameters pRef = base_;
    FDTD1D simRef(pRef);
    simRef.attachMonitor(&mon.incMonitor());
    simRef.run();


    // с диэлектриком
    FDTD1D simWork(base_, layout_);
    simWork.attachMonitor(&mon.totMonitor());
    simWork.run();

    mon.writeTimeSeriesCSV("C:\\Users\\6anna\\CLionProjects\\Sem6Laba3_FDTD\\cmake-build-debug\\monitor_timeseries.csv");

    simWork.writeFieldCSV("C:\\Users\\6anna\\CLionProjects\\Sem6Laba3_FDTD\\cmake-build-debug\\field_Gauss_soft_Silica.csv");


    // Перевод f -> λ:  λ_norm = 1/f,  λ_nm = λ_norm * a_nm
    // Диапазон частот:
    //   fMax = a_nm / lambdaMin  (высокочастотный конец)
    //   fMin = a_nm / lambdaMax  (низкочастотный конец)
    const double fMin = cfgD_.a_nm / cfgD_.lambdaMax;
    const double fMax = cfgD_.a_nm / cfgD_.lambdaMin;
    const int    nF   = cfgD_.nWavelengths;

    auto specF = mon.computeReflection(fMin, fMax, nF);

    // Конвертируем в λ (обращаем порядок: малые f -> большие λ)
    result_.clear();
    result_.reserve(nF);
    for (int i = nF - 1; i >= 0; --i) {
        double f   = specF[i].first;
        double R   = specF[i].second;
        double lam = (f > 1e-30) ? (cfgD_.a_nm / f) : 0.0;
        result_.push_back({lam, R});
    }

    return result_;
}


void DielectricAnalysis::writeCSV(const std::string& filename,
                                   double n1, double n2) const
{
    if (result_.empty()) {
        std::cerr << "DielectricAnalysis::writeCSV: net dannykh, snachala run()\n";
        return;
    }

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "DielectricAnalysis: cannot open " << filename << "\n";
        return;
    }

    const double Rth = fresnelR(n1, n2);

    out << "lambda_nm,R_fdtd,R_theory\n";
    out << std::scientific;
    for (const auto& p : result_)
        out << p.first << "," << p.second << "," << Rth << "\n";

    std::cout << "DielectricAnalysis: written " << result_.size()
              << " points to " << filename << "\n";
    std::cout << "  Theoretical R (Fresnel, n1=" << n1 << ", n2=" << n2
              << ") = " << Rth << "\n";
}