//
// Created by 6anna on 20.04.2026.
//


#include <fstream>
#include <iostream>
#include <iomanip>

#include "FDTD1D.h"

FDTD1D::FDTD1D(const SimulationParameters& p)
    : p_(p),
      Ex_(p.nx + 1, 0.0),
      Hy_(p.nx,     0.0),
      cw_(p.sourceFreq),
      gauss_(p.sourceFreq, p.sourceFWidth)
{}

double FDTD1D::sourceValue(double t) const {
    if (p_.sourceType == SimulationParameters::CW) {
        return cw_(t);
    } else {
        return gauss_(t);
    }
}

void FDTD1D::run() {
    snapshotsEx_.clear();

    // Коэффициенты для однородной немагнитной среды (c = 1)
    const double cE = p_.dt / (p_.eps0 * p_.dx);
    const double cH = p_.dt / (p_.mu0  * p_.dx);

    for (int n = 0; n < p_.numTimeSteps; ++n) {

        // 1) Обновление H^{n+1/2} = H^{n-1/2} - dt/(mu*dx) * (Ex[i+1] - Ex[i])
        for (int i = 0; i < p_.nx; ++i) {
            Hy_[i] -= cH * (Ex_[i + 1] - Ex_[i]);
        }

        // 2) Обновление Ex^{n+1} = Ex^n + dt/(eps*dx) * (Hy[i-1] - Hy[i])
        //    (Знак согласован со стандартной формой ∂E/∂t = -(1/eps) ∂H/∂x)
        for (int i = 1; i < p_.nx; ++i) {
            Ex_[i] -= cE * (Hy_[i] - Hy_[i - 1]);
        }

        // Жёсткие (PEC) границы: Ex[0]=Ex[nx]=0 — будут отражать волну
        Ex_[0]       = 0.0;
        Ex_[p_.nx]   = 0.0;

        // 3) Источник
        const double t = (n + 1) * p_.dt;
        const double s = sourceValue(t);

        if (p_.injectionType == SimulationParameters::SOFT) {
            // Мягкий источник: добавляем значение к Ex в точке
            Ex_[p_.source_pos] += s;
        } else {
            // Через плотность тока J: ∂E/∂t += -J/eps  =>  Ex -= dt/eps * J
            // При J = s получаем инъекцию тока в одну ячейку
            Ex_[p_.source_pos] -= (p_.dt / p_.eps0) * s;
        }

        // Снапшот
        if ((n % p_.snapshotEvery) == 0) {
            snapshotsEx_.push_back(Ex_);
        }
    }

    std::cout << "Simulation finished. Steps: " << p_.numTimeSteps
              << ", snapshots: " << snapshotsEx_.size() << "\n";
}

void FDTD1D::writeFieldCSV(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Cannot open " << filename << " for writing\n";
        return;
    }

    out << std::scientific << std::setprecision(9);
    out << "time_over_fL,x_tilde,Ez\n";

    const std::size_t Nt = snapshotsEx_.size();
    for (std::size_t k = 0; k < Nt; ++k) {
        const double t = (static_cast<double>(k) * p_.snapshotEvery) * p_.dt;
        for (int i = 0; i < p_.nx; ++i) {
            const double x  = (static_cast<double>(i) - p_.source_pos) * p_.dx;
            const double Ez = snapshotsEx_[k][i];
            out << t  << ","
                << x  << ","
                << Ez << "\n";
        }
    }

    out.close();
    std::cout << "Field snapshots written to " << filename << " (CSV)\n";
}