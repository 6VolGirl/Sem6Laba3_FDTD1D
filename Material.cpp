//
// Created by 6anna on 24.04.2026.
//

#include "Material.h"



void MaterialLayout::applyTo(std::vector<double>& epsArr,
         std::vector<double>& muArr) const
{
    for (const auto& r : regions_) {
        for (int i = r.cellStart; i < r.cellEnd && i < static_cast<int>(epsArr.size()); ++i)
            epsArr[i] = r.eps;

        for (int i = r.cellStart; i < r.cellEnd && i < static_cast<int>(muArr.size()); ++i)
            muArr[i] = r.mu;
    }
}


int MaterialLayout::addBraggCavity(int startCell, double nH, double nL, int dH_cells,
                                   int dL_cells, int nPairs, int dCav_cells)
{
    const double epsH = nH * nH;
    const double epsL = nL * nL;
    int cur = startCell;

    // Левое зеркало: (H L)^nPairs
    for (int p = 0; p < nPairs; ++p) {
        regions_.push_back(MaterialRegion(cur, cur + dH_cells, epsH, 1.0,
                                          "H (mirror L)"));
        cur += dH_cells;
        regions_.push_back(MaterialRegion(cur, cur + dL_cells, epsL, 1.0,
                                          "L (mirror L)"));
        cur += dL_cells;
    }

    // Полость (воздух/дефект)
    regions_.push_back(MaterialRegion(cur, cur + dCav_cells, 1.0, 1.0,
                                      "Cavity (air)"));
    int cavStart = cur;
    cur += dCav_cells;

    // Правое зеркало: (H L)^nPairs
    for (int p = 0; p < nPairs; ++p) {
        regions_.push_back(MaterialRegion(cur, cur + dH_cells, epsH, 1.0,
                                          "H (mirror R)"));
        cur += dH_cells;
        regions_.push_back(MaterialRegion(cur, cur + dL_cells, epsL, 1.0,
                                          "L (mirror R)"));
        cur += dL_cells;
    }

    (void)cavStart;
    return cur;
}