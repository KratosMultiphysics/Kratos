// Authors: M.A. Celigueta and S. Latorre (CIMNE)
// Date: October 2015

#include "DEM_D_Linear_viscous_Coulomb_2D_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Linear_viscous_Coulomb2D::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Linear_viscous_Coulomb2D(*this));
        return p_clone;
    }

    void DEM_D_Linear_viscous_Coulomb2D::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_D_Linear_viscous_Coulomb2D to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }


    void DEM_D_Linear_viscous_Coulomb2D::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {

        //Get equivalent Young's Modulus
        const double my_young        = element1->GetYoung();
        const double other_young     = element2->GetYoung();
        const double my_poisson      = element1->GetPoisson();
        const double other_poisson   = element2->GetPoisson();
        const double equiv_young     = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));

        const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        const double other_shear_modulus = 0.5 * other_young / (1.0 + other_poisson);
        const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - other_poisson)/other_shear_modulus);

        //Normal and Tangent elastic constants
        mKn = 2.0 * equiv_young;                                 // 2.0 * equiv_young * sqrt_equiv_radius;
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }

    void DEM_D_Linear_viscous_Coulomb2D::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta) {

        //const double my_radius        = element->GetRadius(); // Get equivalent Radius
        //const double effective_radius = my_radius - ini_delta;
        const double my_young         = element->GetYoung(); // Get equivalent Young's Modulus
        const double walls_young      = wall->GetProperties()[YOUNG_MODULUS];
        const double my_poisson       = element->GetPoisson();
        const double walls_poisson    = wall->GetProperties()[POISSON_RATIO];

        const double equiv_young      = my_young * walls_young / (walls_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - walls_poisson * walls_poisson));

        const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        const double walls_shear_modulus = 0.5 * walls_young / (1.0 + walls_poisson);
        const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - walls_poisson)/walls_shear_modulus);

        /*
        const double effective_young = my_young / (1.0 - my_poisson * my_poisson); // Equivalent Young Modulus for RIGID WALLS!
        const double effective_young = 0.5 * my_young / (1.0 - my_poisson * my_poisson); // Equivalent Young Modulus if the wall has the same E of the sphere
        // Infinite Young was imposed to the wall in the formula that includes both.
        // For flexible walls, the computation is different, Thornton 2012 proposes same Young for both)

        const double effective_shear = 0.5 * my_young / ((1.0 + my_poisson) * (2.0 - my_poisson)); // Equivalent Shear Modulus for RIGID WALLS!
        const double effective_shear = 0.25 * my_young / ((1.0 + my_poisson) * (2.0 - my_poisson)); // Equivalent Shear Modulus if the wall has the same E of the sphere
        // Infinite Shear was imposed to the wall in the formula that includes both.
        // For flexible walls, the computation is different, Thornton 2012 proposes same Shear for both)
        */

        //Normal and Tangent elastic constants
        mKn = 2.0 * equiv_young;                                 // 2.0 * equiv_young * sqrt_equiv_radius;
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }

} /* namespace Kratos.*/
