// System includes
#include <string>
#include <iostream>

// Project includes
//#include "DEM_application.h"
#include "custom_constitutive/DEM_D_JKR_cohesive_law.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEM_D_JKR_Cohesive_Law::DEM_D_JKR_Cohesive_Law() {}

    DEM_D_JKR_Cohesive_Law::~DEM_D_JKR_Cohesive_Law() {}

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_JKR_Cohesive_Law::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_JKR_Cohesive_Law(*this));
        return p_clone;
    }

    std::unique_ptr<DEMDiscontinuumConstitutiveLaw> DEM_D_JKR_Cohesive_Law::CloneUnique() {
        return Kratos::make_unique<DEM_D_JKR_Cohesive_Law>();
    }

    void DEM_D_JKR_Cohesive_Law::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_D_JKR_Cohesive_Law to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    double DEM_D_JKR_Cohesive_Law::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {

        Properties& properties_of_this_contact = element1->GetProperties().GetSubProperties(element2->GetProperties().Id());
        const double equiv_cohesion = properties_of_this_contact[PARTICLE_COHESION];
        const double my_young       = element1->GetYoung();
        const double other_young    = element2->GetYoung();
        const double my_poisson     = element1->GetPoisson();
        const double other_poisson  = element2->GetPoisson();
        const double equiv_young    = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));
        const double my_radius      = element1->GetRadius();
        const double other_radius   = element2->GetRadius();
        const double radius_sum     = my_radius + other_radius;
        const double radius_sum_inv = 1.0 / radius_sum;
        const double equiv_radius   = my_radius * other_radius * radius_sum_inv;
        const double contact_radius = sqrt(equiv_radius * indentation);
        const double cohesive_force = sqrt(8.0 * equiv_cohesion * Globals::Pi * contact_radius * contact_radius * contact_radius * equiv_young);

        return cohesive_force;
    }

    double DEM_D_JKR_Cohesive_Law::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, Condition* const wall, const double indentation) {

        Properties& properties_of_this_contact = element->GetProperties().GetSubProperties(wall->GetProperties().Id());
        const double equiv_cohesion   = properties_of_this_contact[PARTICLE_COHESION];
        const double my_young         = element->GetYoung();
        const double my_poisson       = element->GetPoisson();
        const double equiv_radius     = element->GetRadius(); // Equivalent Radius for RIGID WALLS
        const double walls_young      = wall->GetProperties()[YOUNG_MODULUS];
        const double walls_poisson    = wall->GetProperties()[POISSON_RATIO];
        const double equiv_young      = my_young * walls_young / (walls_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - walls_poisson * walls_poisson));
        //const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        //const double walls_shear_modulus = 0.5 * walls_young / (1.0 + walls_poisson);
        //const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - walls_poisson)/walls_shear_modulus);
        const double contact_radius = sqrt(equiv_radius * indentation);
        const double cohesive_force = sqrt(8.0 * equiv_cohesion * Globals::Pi * contact_radius * contact_radius * contact_radius * equiv_young);

        return cohesive_force;
    }
} // Namespace Kratos
