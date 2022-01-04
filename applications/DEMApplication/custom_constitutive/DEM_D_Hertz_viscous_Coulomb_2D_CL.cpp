// Authors: M.A. Celigueta and S. Latorre (CIMNE)
// Date: July 2015

#include "DEM_D_Hertz_viscous_Coulomb_2D_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Hertz_viscous_Coulomb2D::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Hertz_viscous_Coulomb2D(*this));
        return p_clone;
    }

    std::unique_ptr<DEMDiscontinuumConstitutiveLaw> DEM_D_Hertz_viscous_Coulomb2D::CloneUnique() {
        return Kratos::make_unique<DEM_D_Hertz_viscous_Coulomb2D>();
    }

    void DEM_D_Hertz_viscous_Coulomb2D::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {

        const double my_young      = element1->GetYoung();
        const double other_young   = element2->GetYoung();
        const double my_poisson    = element1->GetPoisson();
        const double other_poisson = element2->GetPoisson();
        const double equiv_young   = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));

        double equiv_poisson;
        if (my_poisson + other_poisson) {
            equiv_poisson = 2.0 * my_poisson * other_poisson / (my_poisson + other_poisson);
        } else {
            equiv_poisson = 0.0;
        }

        mKn = 0.25 * Globals::Pi * equiv_young * 1.0; // This 1.0 is the length (the unitary thickness)
        mKt = mKn * (1.0 - equiv_poisson) / (1.0 - 0.5 * equiv_poisson);
    }

    void DEM_D_Hertz_viscous_Coulomb2D::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta) {

        const double my_young      = element->GetYoung(); //Get equivalent Young's Modulus
        const double walls_young   = wall->GetProperties()[YOUNG_MODULUS];
        const double my_poisson    = element->GetPoisson();
        const double walls_poisson = wall->GetProperties()[POISSON_RATIO];

        const double equiv_young    = my_young * walls_young / (walls_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - walls_poisson * walls_poisson));
        const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        const double walls_shear_modulus = 0.5 * walls_young / (1.0 + walls_poisson);
        const double equiv_shear = 1.0 / ((2.0 - my_poisson) / my_shear_modulus + (2.0 - walls_poisson) / walls_shear_modulus);

        const double element_radius = element->GetRadius();

        mKn = 2.0 * element_radius * equiv_young;
        
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }

    double DEM_D_Hertz_viscous_Coulomb2D::CalculateNormalForce(const double indentation) {
        return mKn * indentation;
    }

} // namespace Kratos
