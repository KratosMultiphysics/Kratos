// System includes
#include <string>
#include <iostream>

// Project includes
#include "../custom_constitutive/DEM_D_DMT_cohesive_law.h"
#include "../custom_elements/spheric_particle.h"

namespace Kratos {

    DEM_D_DMT_Cohesive_Law::DEM_D_DMT_Cohesive_Law() {}

    DEM_D_DMT_Cohesive_Law::~DEM_D_DMT_Cohesive_Law() {}

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_DMT_Cohesive_Law::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_DMT_Cohesive_Law(*this));
        return p_clone;
    }

    void DEM_D_DMT_Cohesive_Law::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_D_DMT_Cohesive_Law to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    double DEM_D_DMT_Cohesive_Law::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {

        const double equiv_cohesion = 0.5 * (element1->GetParticleCohesion() + element2->GetParticleCohesion());
        const double my_radius      = element1->GetRadius();
        const double other_radius   = element2->GetRadius();
        const double radius_sum     = my_radius + other_radius;
        const double radius_sum_inv = 1.0 / radius_sum;
        const double equiv_radius   = my_radius * other_radius * radius_sum_inv;
        const double cohesive_force = 2.0 * Globals::Pi * equiv_cohesion * equiv_radius;

        return cohesive_force;
    }

    double DEM_D_DMT_Cohesive_Law::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, Condition* const wall, const double indentation) {

        const double cohesion       = element->GetParticleCohesion(); // For the time being, this represents the Surface Energy
        const double equiv_cohesion = 0.5 * (cohesion + wall->GetProperties()[WALL_COHESION]);
        const double equiv_radius   = element->GetRadius(); // Equivalent Radius for RIGID WALLS
        const double cohesive_force = 2.0 * Globals::Pi * equiv_cohesion * equiv_radius;

        return cohesive_force;
    }
} // Namespace Kratos
