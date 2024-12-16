/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu
// Date: Feb 2023
/////////////////////////////////////////////////

// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <algorithm>

// Project includes
#include "DEM_parallel_bond_for_membrane_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos{

DEMContinuumConstitutiveLaw::Pointer DEM_parallel_bond_for_membrane::Clone() const{
    DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_parallel_bond_for_membrane(*this));
    return p_clone;
}

//*************************************
// Parameters preparation
//*************************************


// Here we calculate the mKn and mKt
void DEM_parallel_bond_for_membrane::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {

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

        //Get equivalent Shear Modulus
        const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        const double other_shear_modulus = 0.5 * other_young / (1.0 + other_poisson);
        const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - other_poisson)/other_shear_modulus);

        //Literature [Cundall, 2004, "A bonded particle model for rock"] [PFC 7.0 manual]
        mKn = equiv_young * Globals::Pi * min_radius * min_radius / radius_sum;
        mKt = equiv_shear * Globals::Pi * min_radius * min_radius / radius_sum;

}

// TODO: In this function, it is better to replace 'kn_el' with 'kn_bond' and 'kt_el' with 'kt_bond'.
void DEM_parallel_bond_for_membrane::CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                double equiv_poisson, double calculation_area, SphericContinuumParticle* element1, SphericContinuumParticle* element2, double indentation) {

    KRATOS_TRY

    //for bonded part
    const double bond_equiv_young = (*mpProperties)[BOND_YOUNG_MODULUS];
    kn_el = bond_equiv_young * calculation_area / initial_dist;
    kt_el = kn_el / (*mpProperties)[BOND_KNKS_RATIO];

    InitializeContact(element1, element2, indentation);

    KRATOS_CATCH("")

}

//*************************************
// Force calculation
//*************************************

double DEM_parallel_bond_for_membrane::ComputeNormalUnbondedForce(double indentation){
    
    KRATOS_TRY

    if (indentation > 0.0) {
        return mKn * indentation;
    } else {
        return 0.0;
    }

    KRATOS_CATCH("")
}

void DEM_parallel_bond_for_membrane::CalculateUnbondedViscoDampingForce(double LocalRelVel[3],
                                                double UnbondedViscoDampingLocalContactForce[3],
                                                SphericParticle* const element1,
                                                SphericParticle* const element2){
    KRATOS_TRY
    const double my_mass    = element1->GetMass();
    const double other_mass = element2->GetMass();

    const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);

    Properties& properties_of_this_contact = element1->GetProperties().GetSubProperties(element2->GetProperties().Id());
    const double damping_gamma = properties_of_this_contact[DAMPING_GAMMA];

    const double equiv_visco_damp_coeff_normal     = 2.0 * damping_gamma * sqrt(equiv_mass * mKn);
    const double equiv_visco_damp_coeff_tangential = 2.0 * damping_gamma * sqrt(equiv_mass * mKt);

    UnbondedViscoDampingLocalContactForce[0] = - equiv_visco_damp_coeff_tangential * LocalRelVel[0];
    UnbondedViscoDampingLocalContactForce[1] = - equiv_visco_damp_coeff_tangential * LocalRelVel[1];
    UnbondedViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal     * LocalRelVel[2];
    KRATOS_CATCH("")
}


//*************************************
// Moment calculation
//*************************************

void DEM_parallel_bond_for_membrane::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                SphericContinuumParticle* neighbor,
                                                double equiv_young,
                                                double distance,
                                                double calculation_area,
                                                double LocalCoordSystem[3][3],
                                                double ElasticLocalRotationalMoment[3],
                                                double ViscoLocalRotationalMoment[3],
                                                double equiv_poisson,
                                                double indentation) {

    KRATOS_TRY


    ElasticLocalRotationalMoment[0] = 0.0;
    ElasticLocalRotationalMoment[1] = 0.0;
    ElasticLocalRotationalMoment[2] = 0.0;

    ViscoLocalRotationalMoment[0] = 0.0;
    ViscoLocalRotationalMoment[1] = 0.0;
    ViscoLocalRotationalMoment[2] = 0.0;
     
    KRATOS_CATCH("")
}//ComputeParticleRotationalMoments

    
} //namespace Kratos