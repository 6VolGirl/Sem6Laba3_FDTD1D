//
// Created by 6anna on 26.04.2026.
//

#ifndef BRAGGMICROCAVITYANALYSIS_H
#define BRAGGMICROCAVITYANALYSIS_H



#include <vector>
#include <string>
#include <utility>
#include <complex>

#include "SimulationParameters.h"
#include "Material.h"
#include "Monitor.h"
#include "FDTD1D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct CavityResult {
    double lambda_nm;
    double R_cavity;
    double T_cavity;
    double R_periodic;
    double T_periodic;
};

class CavityAnalysis {
public:
    struct ConfigCavity {
        double lambda0_nm = 600.0;

        double nH   = 2.3;
        double nL   = 1.45;
        double nCav = 1.0;

        int nPairs = 5;

        double dx_nm = 1.0;

        double lambdaMin    = 300.0;
        double lambdaMax    = 900.0;
        int    nWavelengths = 400;

        int monRefOffset   = 150;
        int monTransOffset = 50;

        double cavThickness_nm = 0.0;

        std::string fieldPrefix = "cavity_field";
    };

private:
    SimulationParameters      base_;
    ConfigCavity              cfgC_;
    std::vector<CavityResult> result_;

    int structStart_ = 0;
    int structEnd_   = 0;
    int cavStart_    = 0;
    int cavEnd_      = 0;

    void computeLayerCells(int& dH, int& dL, int& dCav) const;

    // Вакуумный прогон: заполняет падающие поля на обоих позициях мониторов
    void runVacuum(int monRefPos, int monTransPos,
                   FieldMonitor& incRef,
                   FieldMonitor& incTrans);

    // Один прогон со структурой
    void runOnce(const MaterialLayout& layout,
                 int monRefPos, int monTransPos,
                 bool saveField,
                 const FieldMonitor& incRef,
                 const FieldMonitor& incTrans,
                 std::vector<double>& outR,
                 std::vector<double>& outT,
                 std::vector<double>& outFreqs);

public:
    explicit CavityAnalysis(const SimulationParameters& base)
        : base_(base), cfgC_() {}

    CavityAnalysis(const SimulationParameters& base, const ConfigCavity& cfgC)
        : base_(base), cfgC_(cfgC) {}

    std::vector<CavityResult> runCavity();

    void writeSpectraCSV(const std::string& filename) const;
    void writeFieldSnapshotsCSV(const std::string& prefix) const;

    const std::vector<CavityResult>& result() const { return result_; }
};






// #include <vector>
// #include <string>
// #include <utility>
//
// #include "SimulationParameters.h"
// #include "Material.h"
// #include "Monitor.h"
// #include "FDTD1D.h"
//
// #ifndef M_PI
// #define M_PI 3.14159265358979323846
// #endif
//
// struct CavityResult {
//     double lambda_nm;
//     double R_cavity;
//     double T_cavity;
//     double R_periodic;  // чисто периодическая (без дефекта)
//     double T_periodic;
// };
//
// class CavityAnalysis {
// public:
//     struct ConfigCavity {
//         double lambda0_nm = 600.0;
//
//         // Материалы зеркал
//         double nH = 2.3;   // TiO2
//         double nL = 1.45;  // SiO2
//         double nCav = 1.0; // полость (воздух)
//
//         int nPairs = 5;
//
//         double dx_nm = 1.0;
//
//         double lambdaMin = 300.0;
//         double lambdaMax = 900.0;
//         int    nWavelengths = 400;
//
//         int monRefOffset   = 150;
//         int monTransOffset = 50;
//
//         // Толщина полости: lambda0/(2*nCav) — полуволновой резонатор
//         // 0 = авторасчёт
//         double cavThickness_nm = 0.0;
//
//         std::string fieldPrefix = "cavity_field";
//     };
//
// private:
//     SimulationParameters base_;
//     ConfigCavity               cfgC_;
//     std::vector<CavityResult> result_;
//
//     std::vector<std::vector<double>> snapsExCav_;
//     std::vector<std::vector<double>> snapsHyCav_;
//     std::vector<int>                 snapsTimes_;
//
//
//     int structStart_ = 0;
//     int structEnd_   = 0;
//     int cavStart_    = 0;
//     int cavEnd_      = 0;
//
//     // Вычисляет толщину слоёв в ячейках из четвертьволнового условия
//     void computeLayerCells(int& dH, int& dL, int& dCav) const;
//
//
//     void runOnce(const MaterialLayout& layout,
//                  int monRefPos, int monTransPos,
//                  bool saveField,
//                  const std::vector<double>& eTransInc,
//                  const std::vector<double>& eHyTransInc,
//                  std::vector<double>& outR,
//                  std::vector<double>& outT,
//                  std::vector<double>& outFreqs);
//
//     void
//     runVacuum(int monRefPos, int monTransPos,
//                    Monitor& mon,
//                    std::vector<double>& eTransIncEx,
//                    std::vector<double>& eTransIncHy);
//
// public:
//     explicit CavityAnalysis(const SimulationParameters& base)
//         : base_(base), cfgC_() {}
//
//     CavityAnalysis(const SimulationParameters& base, const ConfigCavity& cfgC)
//         : base_(base), cfgC_(cfgC) {}
//
//     // Запускает оба прогона (полость + периодическая)
//     std::vector<CavityResult> runCavity();
//
//     // Записывает R/T спектр в CSV
//     void writeSpectraCSV(const std::string& filename) const;
//
//     void writeFieldSnapshotsCSV(const std::string& prefix) const;
//
//     const std::vector<CavityResult>& result() const { return result_; }
//
//
// };



#endif //BRAGGMICROCAVITYANALYSIS_H
