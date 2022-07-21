/////////////////////////////////////////////////
// Author: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: June 2022
/////////////////////////////////////////////////

//TODO: // Here this Quadratic model is established for Parallel Bond Model only.
        // It CAN NOT BE USED AS AN INDEPENDENT CANTACT MODEL UNTIL IT IS COMPLETED.

#include "DEM_D_Quadratic_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Quadratic::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Quadratic(*this));
        return p_clone;
    }

    std::unique_ptr<DEMDiscontinuumConstitutiveLaw> DEM_D_Quadratic::CloneUnique() {
        return Kratos::make_unique<DEM_D_Quadratic>();
    }

    std::string DEM_D_Quadratic::GetTypeOfLaw() {
        std::string type_of_law = "Quadratic";
        return type_of_law;
    }

    void DEM_D_Quadratic::Check(Properties::Pointer pProp) const {

        if(!pProp->Has(ALPHA_K)){
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_ERROR << "Error: Variable ALPHA_K should be present in the properties when using DEM_D_Quadratic_LAW. "<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        //pProp->GetValue(ALPHA_K) = 1.0;
        }
    }


    /////////////////////////
    // DEM-DEM INTERACTION //
    /////////////////////////

    void DEM_D_Quadratic::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {
        
        //Get equivalent Radius
        const double my_radius       = element1->GetRadius();
        const double other_radius    = element2->GetRadius();
        const double radius_sum      = my_radius + other_radius;
        const double radius_sum_inv  = 1.0 / radius_sum;
        const double equiv_radius    = my_radius * other_radius * radius_sum_inv;

        //Get equivalent Young's Modulus
        const double my_young        = element1->GetYoung();
        const double other_young     = element2->GetYoung();
        const double my_poisson      = element1->GetPoisson();
        const double other_poisson   = element2->GetPoisson();
        const double equiv_young     = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));
        const double equiv_poisson   = my_poisson * other_poisson / (my_poisson + other_poisson);

        //Get equivalent Shear Modulus
        const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        const double other_shear_modulus = 0.5 * other_young / (1.0 + other_poisson);
        const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - other_poisson)/other_shear_modulus);

        //Normal and Tangent elastic constants //TODO: UPDATE 
        
        const double alpha_k = element1->GetProperties()[ALPHA_K];
        //const double alpha_k = 4e5; 
        const double aim_radius = std::min(my_radius, other_radius);
        // mKn is from Li,2011 [Modeling of stress-dependent static and dynamic moduli of weak sandstones]
        mKn = alpha_k * equiv_shear * Globals::Pi * aim_radius * aim_radius / radius_sum * (1 - equiv_poisson); 
        mKt = 0.5 * mKn / (1 + equiv_poisson);
    }

    double DEM_D_Quadratic::CalculateNormalForce(const double indentation) {

        return mKn * indentation * indentation;
    }

    void DEM_D_Quadratic::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta) {
        
        //Get effective Radius
        const double my_radius           = element->GetRadius(); //Get equivalent Radius
        const double effective_radius    = my_radius - ini_delta;

        //Get equivalent Young's Modulus
        const double my_young            = element->GetYoung();
        const double walls_young         = wall->GetProperties()[YOUNG_MODULUS];
        const double my_poisson          = element->GetPoisson();
        const double walls_poisson       = wall->GetProperties()[POISSON_RATIO];
        const double equiv_young         = my_young * walls_young / (walls_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - walls_poisson * walls_poisson));
        const double equiv_poisson       = my_poisson * walls_poisson / (my_poisson + walls_poisson);

        //Get equivalent Shear Modulus
        const double my_shear_modulus    = 0.5 * my_young / (1.0 + my_poisson);
        const double walls_shear_modulus = 0.5 * walls_young / (1.0 + walls_poisson);
        const double equiv_shear         = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - walls_poisson)/walls_shear_modulus);
 
        mKn = equiv_young * Globals::Pi * my_radius * my_radius / (my_radius + my_radius); 
        mKt = equiv_shear * Globals::Pi * my_radius * my_radius / (my_radius + my_radius);
    }

} //namespace Kratos
