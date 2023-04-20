/////////////////////////////////////////////////
// Author: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: June 2022
/////////////////////////////////////////////////

#include "DEM_D_Linear_classic_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos{

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Linear_classic::Clone() const {
        
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Linear_classic(*this));
        return p_clone;
    }

    std::unique_ptr<DEMDiscontinuumConstitutiveLaw> DEM_D_Linear_classic::CloneUnique() {
        return Kratos::make_unique<DEM_D_Linear_classic>();
    }

    void DEM_D_Linear_classic::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {

        //Get equivalent Radius
        const double my_radius       = element1->GetRadius();
        const double other_radius    = element2->GetRadius();
        const double radius_sum      = my_radius + other_radius;
        const double min_radius      = std::min(my_radius, other_radius);

        //Get equivalent Young's Modulus
        const double my_young        = element1->GetYoung();
        const double other_young     = element2->GetYoung();
        const double my_poisson      = element1->GetPoisson();
        const double other_poisson   = element2->GetPoisson();
        const double equiv_young     = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));
        const double equiv_poisson   = 2.0 * my_poisson * other_poisson / (my_poisson + other_poisson);

        //Get equivalent Shear Modulus
        //const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        //const double other_shear_modulus = 0.5 * other_young / (1.0 + other_poisson);
        //const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - other_poisson)/other_shear_modulus);

        //Literature [Cundall, 2004, "A bonded particle model for rock"] [PFC 7.0 manual]
        mKn = equiv_young * Globals::Pi * min_radius * min_radius / radius_sum;
        mKt = mKn / (2.0 * (equiv_poisson + 1.0));
    }

    void DEM_D_Linear_classic::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta) {
        
        //Get effective Radius
        const double my_radius           = element->GetRadius(); //Get equivalent Radius
        const double effective_radius    = my_radius - ini_delta;

        //Get equivalent Young's Modulus
        const double my_young            = element->GetYoung();
        const double walls_young         = wall->GetProperties()[YOUNG_MODULUS];
        const double my_poisson          = element->GetPoisson();
        const double walls_poisson       = wall->GetProperties()[POISSON_RATIO];
        const double equiv_young         = my_young * walls_young / (walls_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - walls_poisson * walls_poisson));
        const double equiv_poisson       = 2.0 * my_poisson * walls_poisson / (my_poisson + walls_poisson);

        //Get equivalent Shear Modulus
        //const double my_shear_modulus    = 0.5 * my_young / (1.0 + my_poisson);
        //const double walls_shear_modulus = 0.5 * walls_young / (1.0 + walls_poisson);
        //const double equiv_shear         = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - walls_poisson)/walls_shear_modulus);

        //Literature [Cundall, 2004, "A bonded particle model for rock"]
        mKn = equiv_young * Globals::Pi * effective_radius;
        mKt = mKn / (2.0 * (equiv_poisson + 1.0));
    }
} // namespace Kratos