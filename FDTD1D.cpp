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
      eps_(p.nx + 1, p.eps0),
      mu_(p.nx,     p.mu0),
      pml_(p.nx, p.pmlThickness, p.eps0, p.mu0, p.pmlDamping, p.pmlProfilePower, p.dx),
      cw_(p.sourceFreq),
      gauss_(p.sourceFreq, p.sourceFWidth)
{
    sigmaE_ = pml_.sigmaE;
    sigmaM_ = pml_.sigmaM;

}

FDTD1D::FDTD1D(const SimulationParameters& p, const MaterialLayout& layout)
    : p_(p),
      Ex_(p.nx + 1, 0.0),
      Hy_(p.nx,     0.0),
      eps_(p.nx + 1, p.eps0),   // сначала заполняем вакуумом
      mu_ (p.nx,     p.mu0),
      pml_(p.nx, p.pmlThickness, p.eps0, p.mu0,
           p.pmlDamping, p.pmlProfilePower, p.dx),
      cw_(p.sourceFreq),
      gauss_(p.sourceFreq, p.sourceFWidth)
{
    sigmaE_ = pml_.sigmaE;
    sigmaM_ = pml_.sigmaM;

    layout.applyTo(eps_, mu_);
}

double FDTD1D::sourceValue(double t) const {
    if (p_.sourceType == SimulationParameters::CW) {
        return cw_(t);
    } else {
        return gauss_(t);
    }
}

void FDTD1D::run() {
    snapshotsEx_.clear();

    if (monitor_) monitor_->clear();

    //Учёт sigmaM
    std::vector<double> Da(p_.nx), Db(p_.nx);
    for (int i = 0; i < p_.nx; ++i) {
        const double denom = 2.0 * mu_[i] + sigmaM_[i] * p_.dt;
        Da[i] = (2.0 * mu_[i] - sigmaM_[i] * p_.dt) / denom;
        Db[i] = (2.0 * p_.dt) / (p_.dx * denom);
    }

    //Учёт sigmaE
    std::vector<double> Ca(p_.nx + 1), Cb(p_.nx + 1);
    for (int i = 0; i <= p_.nx; ++i) {
        const double denom = 2.0 * eps_[i] + sigmaE_[i] * p_.dt;
        Ca[i] = (2.0 * eps_[i] - sigmaE_[i] * p_.dt) / denom;
        Cb[i] = (2.0 * p_.dt) / (p_.dx * denom);
    }

    for (int n = 0; n < p_.numTimeSteps; ++n) {

        // 1) Обновление H^{n+1/2} = H^{n-1/2} - dt/(mu*dx) * (Ex[i+1] - Ex[i])
        for (int i = 0; i < p_.nx; ++i) {
            Hy_[i] = Da[i] * Hy_[i] - Db[i] * (Ex_[i + 1] - Ex_[i]);
        }

        // 2) Обновление Ex^{n+1} = Ex^n + dt/(eps*dx) * (Hy[i-1] - Hy[i])
        for (int i = 1; i < p_.nx; ++i) {
            Ex_[i] = Ca[i] * Ex_[i] - Cb[i] * (Hy_[i] - Hy_[i - 1]);
        }

        // Границы: Ex[0]=Ex[nx]=0
        Ex_[0]       = 0.0;
        Ex_[p_.nx]   = 0.0;

        // 3) Источник
        const double t = (n + 1) * p_.dt;
        const double s = sourceValue(t);

        if (p_.injectionType == SimulationParameters::SOFT) {
            // Мягкий источник
            Ex_[p_.source_pos] += s;
        } else {
            // Через плотность тока J: ∂E/∂t += -J/eps
            Ex_[p_.source_pos] -= (p_.dt / p_.eps0) * s;
        }

        if (monitor_) {
            monitor_->record(Ex_);
        }

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
