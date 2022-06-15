/////////////////////////////////////////////////
// Author: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: June 2022
/////////////////////////////////////////////////

// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <algorithm>

// Project includes
#include "DEM_parallel_bond_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos{

//TODO: understand this
DEMContinuumConstitutiveLaw::Pointer DEM_parallel_bond::Clone() const{
    DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_parallel_bond(*this));
    return p_clone;
}

//*************************************
// Parameters preparation
//*************************************

void DEM_parallel_bond::TransferParametersToProperties(const Parameters& parameters, Properties::Pointer pProp){
    BaseClassType::TransferParametersToProperties(parameters, pProp);

    //TODO: add the parameters need to be transferred

}

void DEM_parallel_bond::Check(Properties::Pointer pProp) const {
    
    //How many parameters do I need?
    //two parts: discontinuum part and continuum part

    //*********discontinuum part **************
    //for [DEM_D_Hertz_viscous_Coulomb]
    //here we only mention static_friction and dynamic_friction 

    if(!pProp->Has(STATIC_FRICTION)){
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable STATIC_FRICTION should be present in the properties when using DEMContinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(STATIC_FRICTION) = 0.0;
    }

    if(!pProp->Has(DYNAMIC_FRICTION)){
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable DYNAMIC_FRICTION should be present in the properties when using DEMContinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(DYNAMIC_FRICTION) = 0.0;
    }

    if(!pProp->Has(FRICTION_DECAY)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable FRICTION_DECAY should be present in the properties when using DEMContinuumConstitutiveLaw. 500.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(FRICTION_DECAY) = 500.0;
    }

    if(!pProp->Has(COEFFICIENT_OF_RESTITUTION)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable COEFFICIENT_OF_RESTITUTION should be present in the properties when using DEMContinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(COEFFICIENT_OF_RESTITUTION) = 0.0;
    }

    if(!pProp->Has(ROLLING_FRICTION)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable ROLLING_FRICTION should be present in the properties when using DEMContinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(ROLLING_FRICTION) = 0.0;
    }

    if(!pProp->Has(ROLLING_FRICTION_WITH_WALLS)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable ROLLING_FRICTION_WITH_WALLS should be present in the properties when using DEMContinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(ROLLING_FRICTION_WITH_WALLS) = 0.0;
    }

    //********** continuum part ***************
    //for [DEM_parallel_bond_CL]

    if(!pProp->Has(BOND_YOUNG_MODULUS)) {

        KRATOS_ERROR << "Variable BONDED_MATERIAL_YOUNG_MODULUS should be present in the properties when using DEM_parallel_bond_CL."<<std::endl;

    }

    if(!pProp->Has(BOND_KNKS_RATIO)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable BOND_KNKS_RATIO should be present in the properties when using DEM_parallel_bond_CL. 2.5 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(BOND_KNKS_RATIO) = 2.5;
    }

    if(!pProp->Has(BOND_SIGMA_MAX)) {

        KRATOS_ERROR << "Variable BOND_SIGMA_MAX should be present in the properties when using DEM_parallel_bond_CL."<<std::endl;

    }

    if(!pProp->Has(BOND_SIGMA_MAX_DEVIATION)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable BOND_SIGMA_MAX_DEVIATION should be present in the properties when using DEM_parallel_bond_CL. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(BOND_SIGMA_MAX_DEVIATION) = 0.0;
    }

    if(!pProp->Has(BOND_TAU_ZERO)) {

        KRATOS_ERROR << "Variable BOND_TAU_ZERO should be present in the properties when using DEM_parallel_bond_CL."<<std::endl;

    }

    if(!pProp->Has(BOND_TAU_ZERO_DEVIATION)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable BOND_TAU_ZERO_DEVIATION should be present in the properties when using DEM_parallel_bond_CL. 0.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(BOND_TAU_ZERO_DEVIATION) = 0.0;
    }

    if(!pProp->Has(BOND_INTERNAL_FRICC)) {

        KRATOS_ERROR << "Variable BOND_INTERNAL_FRICC should be present in the properties when using DEM_parallel_bond_CL."<<std::endl;

    }

    if(!pProp->Has(BOND_ROTATIONAL_MOMENT_COEFFICIENT)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable BOND_ROTATIONAL_MOMENT_COEFFICIENT should be present in the properties when using DEM_parallel_bond_CL. 0.1 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(BOND_ROTATIONAL_MOMENT_COEFFICIENT) = 0.1;
    }

    if(!pProp->Has(BOND_RADIUS_FACTOR)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable BOND_RADIUS_FACTOR should be present in the properties when using DEM_parallel_bond_CL. 1.0 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(BOND_RADIUS_FACTOR) = 1.0;
    }

    if (!pProp->Has(IS_UNBREAKABLE)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable IS_UNBREAKABLE was not present in the properties when using DEM_parallel_bond_CL. False value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(IS_UNBREAKABLE) = false;
    }
} // CHECK()

// TODO: Calculate Contact Area = Calculate Bond area?
void DEM_parallel_bond::CalculateContactArea(double radius, double other_radius, double& calculation_area) {

    KRATOS_TRY

    double aim_radius = min(radius, other_radius);
    calculation_area = Globals::Pi * aim_radius * aim_radius;

    KRATOS_CATCH("")
}

double DEM_parallel_bond::CalculateContactArea(double radius, double other_radius, Vector& v) {
    
    KRATOS_TRY

    double a = 0.0;
    CalculateContactArea(radius, other_radius, a);
    unsigned int old_size = v.size();
    Vector backup = v;
    v.resize(old_size + 1, false);
    v[old_size] = a;
    for (unsigned int i=0; i<old_size; i++) {
        v[i] = backup[i];
    }
    return a;

    KRATOS_CATCH("")
}

void DEM_parallel_bond::GetcontactArea(const double radius, const double other_radius, const Vecotr& vector_of_initial_areas, const int neighbour_position, double& calculation_area) {
    
    KRATOS_TRY

    if(vector_of_initial_areas.size()){
        calculation_area = vector_of_initial_areas[neighbour_position];
    }
    else{
        CalculateContactArea(radius, other_radius, calculation_area);
    }

    KRATOS_CATCH("")
}

// TODO: In this function, it is better to replace 'kn_el' with 'kn_bond' and 'kt_el' with 'kt_bond'.
void DEM_parallel_bond::CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                double equiv_poisson, double calculation_area, SphericContinuumParticle* element1, SphericContinuumParticle* element2, double indentation) {

    KRATOS_TRY

    //for unbonded part
    //TODO: check whether the unbonded elastic parameters should be calculated here
    //maybe they should not be here.

    //for bonded part
    const double bond_equiv_young = (*mpProperties)[BOND_YOUNG_MODULUS];
    kn_el = bond_equiv_young * calculation_area / initial_dist;
    kt_el = kn_el / (*mpProperties)[BOND_KNKS_RATIO];

    KRATOS_CATCH("")

}

double DEM_parallel_bond::LocalMaxSearchDistance(const int i,
                            SphericContinuumParticle* element1,
                            SphericContinuumParticle* element2) {
    KRATOS_TRY

    // TODO: maybe this function is unnecessary

    KRATOS_CATCH("")
}

double DEM_parallel_bond::GetContactSigmaMax(){

    KRATOS_TRY

    // TODO: maybe this function is unnecessary

    KRATOS_CATCH("")    
}



//*************************************
// Bonded force calculation
//*************************************
void DEM_parallel_bond::CalculateForces(const ProcessInfo& r_process_info,
                            double OldLocalElasticContactForce[3],
                            double LocalElasticContactForce[3],
                            double LocalElasticExtraContactForce[3],
                            double LocalCoordSystem[3][3],
                            double LocalDeltDisp[3],
                            const double kn_el,
                            const double kt_el,
                            double& contact_sigma,
                            double& contact_tau,
                            double& failure_criterion_state,
                            double equiv_young,
                            double equiv_shear,
                            double indentation,
                            double calculation_area,
                            double& acumulated_damage,
                            SphericContinuumParticle* element1,
                            SphericContinuumParticle* element2,
                            int i_neighbour_count,
                            int time_steps,
                            bool& sliding,
                            double &equiv_visco_damp_coeff_normal,
                            double &equiv_visco_damp_coeff_tangential,
                            double LocalRelVel[3],
                            double ViscoDampingLocalContactForce[3]) {
                           
    KRATOS_TRY

    CalculateNormalForces(LocalElasticContactForce,
                        kn_el,
                        equiv_young,
                        indentation,
                        calculation_area,
                        acumulated_damage,
                        element1,
                        element2,
                        i_neighbour_count,
                        time_steps,
                        r_process_info);
            
    CalculateViscoDampingCoeff(equiv_visco_damp_coeff_normal,
                            equiv_visco_damp_coeff_tangential,
                            element1,
                            element2,
                            kn_el,
                            kt_el);

    CalculateViscoDamping(LocalRelVel,
                        ViscoDampingLocalContactForce,
                        indentation,
                        equiv_visco_damp_coeff_normal,
                        equiv_visco_damp_coeff_tangential,
                        sliding,
                        element1->mIniNeighbourFailureId[i_neighbour_count]);

    // Tangential forces are calculated after the viscodamping because the frictional limit bounds the sum of elastic plus viscous forces
    CalculateTangentialForces(OldLocalElasticContactForce,
                            LocalElasticContactForce,
                            LocalElasticExtraContactForce,
                            ViscoDampingLocalContactForce,
                            LocalCoordSystem,
                            LocalDeltDisp,
                            LocalRelVel,
                            kt_el,
                            equiv_shear,
                            contact_sigma,
                            contact_tau,
                            indentation,
                            calculation_area,
                            failure_criterion_state,
                            element1,
                            element2,
                            i_neighbour_count,
                            sliding,
                            r_process_info);

    KRATOS_CATCH("") 
}

void DEM_parallel_bond::CalculateNormalForces(double LocalElasticContactForce[3],
                const double kn_el,
                double equiv_young,
                double indentation,
                double calculation_area,
                double& acumulated_damage,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                int i_neighbour_count,
                int time_steps,
            const ProcessInfo& r_process_info) {

    KRATOS_TRY

    int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
    double BondedLocalElasticContactForce2 = 0.0;
    const double bonded_indentation = indentation - mInitialIndentationForBondedPart;                                                                                                          
    double unbonded_indentation = indentation - element1->GetInitialDelta(i_neighbour_count);
    const double& bond_sigma_max = (*mpProperties)[BOND_SIGMA_MAX];
    double bond_current_sigma = 0.0;
    const double& bond_rotational_moment_coefficient =(*mpProperties)[BOND_ROTATIONAL_MOMENT_COEFFICIENT];

    if (!failure_type){ //if the bond is not broken
        BondedLocalElasticContactForce2 = kn_el * bonded_indentation;
    } else { //else the bond is broken
        BondedLocalElasticContactForce2 = 0.0;
    }

    
    ComputeNormalUnbondedForce(unbonded_indentation);

    KRATOS_CATCH("")  

} // CalculateNormalForces


//Moment calculation


//Bond Failure check



//*************************************
// unbonded froce calculation
//*************************************


    
} //namespace Kratos