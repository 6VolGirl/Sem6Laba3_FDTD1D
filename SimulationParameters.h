//
// Created by 6anna on 20.04.2026.
//

#ifndef SIMULATIONPARAMETERS_H
#define SIMULATIONPARAMETERS_H


struct SimulationParameters {
    int    nx            = 1000;   // число ячеек по Ex
    int    numTimeSteps  = 2000;   // число временных шагов
    double dx            = 1.0;    // шаг по пространству (норм)
    double dt            = 0.5;    // шаг по времени (c=1)
    double courantNumber = 0.5;

    // Параметры среды
    double eps0 = 1.0;
    double mu0  = 1.0;

    // Источник
    int    source_pos    = 200;
    double sourceFreq    = 0.05;   // частота
    double sourceFWidth  = 0.02;   // ширина спектра для гауссова импульса

    // Тип источника
    enum SourceType  { CW = 0, GAUSS = 1 };
    enum InjectionType { SOFT = 0, CURRENT = 1 };

    SourceType   sourceType    = GAUSS;
    InjectionType injectionType = SOFT;

    // PML
    int    pmlThickness   = 20;
    double pmlDamping     = 1e-8;
    int    pmlProfilePower = 3;

    int snapshotEvery = 2;
};



#endif //SIMULATIONPARAMETERS_H
