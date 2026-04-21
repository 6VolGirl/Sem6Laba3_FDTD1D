//
// Created by 6anna on 21.04.2026.
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