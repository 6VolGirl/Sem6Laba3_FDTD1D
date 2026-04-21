//
// Created by 6anna on 20.04.2026.
//

#ifndef PMLCOEFFICIENTS_H
#define PMLCOEFFICIENTS_H




#include <vector>

class PMLSigma {
public:

    std::vector<double> sigmaE;
    std::vector<double> sigmaM;


    // pmlCells  - толщина одного PML-слоя в ячейках
    // damping   - желаемое остаточное отражение амплитуды delta << 1
    // power     - степень профиля m: 0 (constant), 1 (linear), 2, 3, ...
    PMLSigma(int nx, int pmlCells, double eps, double mu, double damping,
             int power, double dx);
};



#endif //PMLCOEFFICIENTS_H
