//
// Created by 6anna on 20.04.2026.
//

#include <cmath>
#include "PMLCoefficients.h"

PMLSigma::PMLSigma(int nx, int pmlCells, double eps, double mu,
                   double damping, int power, double dx)
    : sigmaE(nx + 1, 0.0), sigmaM(nx, 0.0)
{
    if (pmlCells <= 0 || damping >= 1.0 || damping <= 0.0) {
        return;
    }

    const double eta = std::sqrt(mu / eps);        // волновое сопротивление
    const double L   = pmlCells * dx;              // физическая толщина PML
    const int    m   = power;                      // степень профиля

    const double sigmaMax = -(m + 1.0) * std::log(damping) / (eta * L);

    // profile(u) : u = расстояние вглубь PML, нормированное на толщину, 0..1
    auto profile = [&](double u) -> double {
        if (u <= 0.0) return 0.0;
        if (u >= 1.0) u = 1.0;
        if (m == 0) return sigmaMax;
        return sigmaMax * std::pow(u, m);
    };


    for (int i = 0; i <= nx; ++i) {
        if (i < pmlCells) {
            double u = double(pmlCells - i) / double(pmlCells);
            sigmaE[i] = profile(u);
        } else if (i > nx - pmlCells) {
            double u = double(i - (nx - pmlCells)) / double(pmlCells);
            sigmaE[i] = profile(u);
        }
    }

    for (int i = 0; i < nx; ++i) {
        double xh = i + 0.5;
        double sig = 0.0;
        if (xh < pmlCells) {
            double u = (pmlCells - xh) / double(pmlCells);
            sig = profile(u);
        } else if (xh > nx - pmlCells) {
            double u = (xh - (nx - pmlCells)) / double(pmlCells);
            sig = profile(u);
        }
        sigmaM[i] = sig * mu / eps;
    }
}