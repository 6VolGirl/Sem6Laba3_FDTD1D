//
// Created by 6anna on 20.04.2026.
//

#ifndef SIMULATIONPARAMETERS_H
#define SIMULATIONPARAMETERS_H


struct SimulationParameters {
    // Сетка
    int    nx            = 1000;
    int    numTimeSteps  = 2000;
    double dx            = 1.0;
    double dt            = 0.5;
    double courantNumber = 0.5;

    //(c=1)
    double eps0 = 1.0;
    double mu0  = 1.0;

    // Источник
    int    source_pos    = 200;
    double sourceFreq    = 2.0;
    double sourceFWidth  = 1.9;

    enum SourceType    { CW = 0, GAUSS = 1 };
    enum InjectionType { SOFT = 0, CURRENT = 1 };

    SourceType    sourceType    = GAUSS;
    InjectionType injectionType = SOFT;

    // PML
    int    pmlThickness   = 20;      // толщина одного слоя в ячейках
    double pmlDamping     = 1e-8;    // delta — желаемое остаточное отражение амплитуды
    int    pmlProfilePower = 3;      // m: 0 (const), 1 (linear), 2, 3, ...

    int snapshotEvery = 2;
};



#endif //SIMULATIONPARAMETERS_H
