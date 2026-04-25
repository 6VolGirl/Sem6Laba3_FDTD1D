//
// Created by 6anna on 24.04.2026.
//

#ifndef MATERIAL_H
#define MATERIAL_H


#include <vector>
#include <string>


struct MaterialRegion {
    int    cellStart;
    int    cellEnd;
    double eps;
    double mu;

    std::string name;


    MaterialRegion(int start, int end, double eps_, double mu_ = 1.0,
                   std::string nm = "")
        : cellStart(start), cellEnd(end), eps(eps_), mu(mu_),
          name(std::move(nm)) {}


    // Кварц: n = 1.45, eps = n^2 = 2.1025, mu = 1
    static MaterialRegion Silica(int start, int end) {
        return MaterialRegion(start, end, 1.45 * 1.45, 1.0, "SiO2 (n=1.45)");
    }

    static MaterialRegion Vacuum(int start, int end) {
        return MaterialRegion(start, end, 1.0, 1.0, "Vacuum");
    }
};


//  Коллекция материальных областей
class MaterialLayout {
private:
    std::vector<MaterialRegion> regions_;

public:
    void add(const MaterialRegion& r) { regions_.push_back(r); }

    // Заполнить массивы eps[nx+1] и mu[nx]
    void applyTo(std::vector<double>& epsArr,
                 std::vector<double>& muArr) const;

    const std::vector<MaterialRegion>& regions() const { return regions_; }

    bool empty() const { return regions_.empty(); }


};



#endif //MATERIAL_H
