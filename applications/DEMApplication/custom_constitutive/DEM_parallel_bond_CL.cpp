/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: June 2022
/////////////////////////////////////////////////

// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <algorithm>
//#include <iomanip>

// Project includes
#include "dem_contact.h"
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

void DEM_parallel_bond::Initialize(SphericContinuumParticle* element1,
                                   SphericContinuumParticle* element2,
                                   Properties::Pointer pProps) {
    mpProperties = pProps;
    mDebugPrintingOption = false;
    if (!mpProperties->Has(DEBUG_PRINTING_OPTION)) {
        mDebugPrintingOption = false;
    } else {
        mDebugPrintingOption = bool((*mpProperties)[DEBUG_PRINTING_OPTION]);
    }
    if (mDebugPrintingOption) {
        if (!mpProperties->Has(DEBUG_PRINTING_ID_1) || !mpProperties->Has(DEBUG_PRINTING_ID_2)) {
            KRATOS_WARNING("DEM") << "\nWARNING: We are currently in DEBUG PRINTING mode, so the ids of the two particles involved must be given.\n\n";
        }
    }
}

void DEM_parallel_bond::TransferParametersToProperties(const Parameters& parameters, Properties::Pointer pProp){
    BaseClassType::TransferParametersToProperties(parameters, pProp);
}

std::string DEM_parallel_bond::GetTypeOfLaw() {
        std::string type_of_law = "parallel_bond_CL";
        return type_of_law;
    }

void DEM_parallel_bond::Check(Properties::Pointer pProp) const {
    
    //two parts: discontinuum part and continuum part

    //*********discontinuum part **************

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

    if(!pProp->Has(BOND_ROTATIONAL_MOMENT_COEFFICIENT_NORMAL)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable BOND_ROTATIONAL_MOMENT_COEFFICIENT_NORMAL should be present in the properties when using DEM_parallel_bond_CL. 0.1 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(BOND_ROTATIONAL_MOMENT_COEFFICIENT_NORMAL) = 0.1;
    }

    if(!pProp->Has(BOND_ROTATIONAL_MOMENT_COEFFICIENT_TANGENTIAL)) {
        KRATOS_WARNING("DEM")<<std::endl;
        KRATOS_WARNING("DEM")<<"WARNING: Variable BOND_ROTATIONAL_MOMENT_COEFFICIENT_TANGENTIAL should be present in the properties when using DEM_parallel_bond_CL. 0.1 value assigned by default."<<std::endl;
        KRATOS_WARNING("DEM")<<std::endl;
        pProp->GetValue(BOND_ROTATIONAL_MOMENT_COEFFICIENT_TANGENTIAL) = 0.1;
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

void DEM_parallel_bond::CalculateContactArea(double radius, double other_radius, double& calculation_area) {

    KRATOS_TRY

    const double bond_radius_factor = (*mpProperties)[BOND_RADIUS_FACTOR];
    double bond_radius = std::min(radius, other_radius) * bond_radius_factor;
    calculation_area = Globals::Pi * bond_radius * bond_radius;

    KRATOS_CATCH("")
}

//TODO: it seems this function never been called?
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

void DEM_parallel_bond::GetContactArea(const double radius, const double other_radius, const Vector& vector_of_initial_areas, const int neighbour_position, double& calculation_area) {
    
    KRATOS_TRY

    CalculateContactArea(radius, other_radius, calculation_area);

    KRATOS_CATCH("")
}

// Here we calculate the mKn and mKt
void DEM_parallel_bond::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {

}

// TODO: In this function, it is better to replace 'kn_el' with 'kn_bond' and 'kt_el' with 'kt_bond'.
void DEM_parallel_bond::CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                double equiv_poisson, double calculation_area, SphericContinuumParticle* element1, SphericContinuumParticle* element2, double indentation) {

    KRATOS_TRY

    //for bonded part
    const double bond_equiv_young = (*mpProperties)[BOND_YOUNG_MODULUS];
    kn_el = bond_equiv_young * calculation_area / initial_dist;
    kt_el = kn_el / (*mpProperties)[BOND_KNKS_RATIO];

    InitializeContact(element1, element2, indentation);

    KRATOS_CATCH("")

}

double DEM_parallel_bond::LocalMaxSearchDistance(const int i,
                            SphericContinuumParticle* element1,
                            SphericContinuumParticle* element2) {
    KRATOS_TRY

    double tension_limit;

    // calculation of equivalent Young modulus
    const double& equiv_young = (*mpProperties)[BOND_YOUNG_MODULUS];

    const double my_radius = element1->GetRadius();
    const double other_radius = element2->GetRadius();
    double calculation_area = 0.0;

    Vector& vector_of_contact_areas = element1->GetValue(NEIGHBOURS_CONTACT_AREAS);
    GetContactArea(my_radius, other_radius, vector_of_contact_areas, i, calculation_area);

    double radius_sum = my_radius + other_radius;
    double initial_delta = element1->GetInitialDelta(i);
    double initial_dist = radius_sum - initial_delta;

    // calculation of elastic constants
    double kn_el = equiv_young * calculation_area / initial_dist;

    //tension_limit = GetContactSigmaMax();
    tension_limit = (*mpProperties)[BOND_SIGMA_MAX]; //TODO: add BOND_SIGMA_MAX_DEVIATION

    const double max_normal_bond_force = tension_limit * calculation_area;
    double u_max = max_normal_bond_force / kn_el;

    //TODO: need to choose whether the [if] below is necessary
    if (u_max > 2.0 * radius_sum) {
        u_max = 2.0 * radius_sum;
    } // avoid error in special cases with too high tensile
    

    return u_max;

    KRATOS_CATCH("")
}

//*************************************
// Force calculation
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
                        r_process_info,
                        contact_sigma);
            
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
                        element1->mIniNeighbourFailureId[i_neighbour_count],
                        i_neighbour_count,
                        element1,
                        element2);

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


double DEM_parallel_bond::ComputeNormalUnbondedForce(double indentation){
    
    KRATOS_TRY

    KRATOS_ERROR << "This function shouldn't be accessed here, use basic contact model instead."<<std::endl
                << "Maybe you are using \"DEM_parallel_bond\" for DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME." <<std::endl
                << "Unfortunately, you can only input one of the names listed below." <<std::endl
                << "1. DEM_parallel_bond_Linear" <<std::endl 
                << "2. DEM_parallel_bond_Hertz" <<std::endl
                << "3. DEM_parallel_bond_Quadratic" <<std::endl;
    /*
    if (unbonded_indentation > 0.0) {
        return mKn * unbonded_indentation;
    } else {
        return 0.0;
    }
    */

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
                const ProcessInfo& r_process_info,
                double& contact_sigma) {

    KRATOS_TRY
    
    int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
    const double bonded_indentation = indentation - mInitialIndentationForBondedPart;                                                                                                          

    double BondedLocalElasticContactForce2 = 0.0;

    if (!failure_type){ //if the bond is not broken
        BondedLocalElasticContactForce2 = kn_el * bonded_indentation;
    } else { //else the bond is broken
        //if the bond is broken, we still calculate the normal compressive force but not the normal tensile force
        if (bonded_indentation > 0.0){
            BondedLocalElasticContactForce2 = kn_el * bonded_indentation;
        } else {
            BondedLocalElasticContactForce2 = 0.0;
        }
    }

    if (indentation > 0.0) {
        mUnbondedLocalElasticContactForce2 = ComputeNormalUnbondedForce(indentation);
    } else {
        mUnbondedLocalElasticContactForce2 = 0.0;
    }
        
    if(calculation_area){
        contact_sigma = BondedLocalElasticContactForce2 / calculation_area;
    }

    LocalElasticContactForce[2] = BondedLocalElasticContactForce2 + mUnbondedLocalElasticContactForce2;

    if (LocalElasticContactForce[2]) {
        mBondedScalingFactor[2] = BondedLocalElasticContactForce2 / LocalElasticContactForce[2]; 
    } else {
        mBondedScalingFactor[2] = 0.0;
    }

    //for debug
    if (mDebugPrintingOption) {

        const long unsigned int& sphere_id = (*mpProperties)[DEBUG_PRINTING_ID_1];
        const long unsigned int& neigh_sphere_id = (*mpProperties)[DEBUG_PRINTING_ID_2];

        if ((element1->Id() == sphere_id) && (element2->Id() == neigh_sphere_id)) {
            std::ofstream normal_forces_file("delta_stress_normal.txt", std::ios_base::out | std::ios_base::app);
            normal_forces_file << r_process_info[TIME] << " " << bonded_indentation/*2*/ << " " << contact_sigma/*3*/ << " " << indentation <<'\n'; 
            normal_forces_file.flush();
            normal_forces_file.close();
        }
    }

    KRATOS_CATCH("")  

} // CalculateNormalForces

void DEM_parallel_bond::CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                        double &equiv_visco_damp_coeff_tangential,
                        SphericContinuumParticle* element1,
                        SphericContinuumParticle* element2,
                        const double kn_el,
                        const double kt_el) {

    KRATOS_TRY

    const double my_mass    = element1->GetMass();
    const double other_mass = element2->GetMass();

    const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);

    // TODO: check the value of gamma
    const double equiv_gamma = (*mpProperties)[DAMPING_GAMMA];

    equiv_visco_damp_coeff_normal     = 2.0 * equiv_gamma * sqrt(equiv_mass * kn_el);
    equiv_visco_damp_coeff_tangential = 2.0 * equiv_gamma * sqrt(equiv_mass * kt_el);

    KRATOS_CATCH("")                
}

void DEM_parallel_bond::CalculateUnbondedViscoDampingForce(double LocalRelVel[3],
                                                double UnbondedViscoDampingLocalContactForce[3],
                                                SphericParticle* const element1,
                                                SphericParticle* const element2){
    KRATOS_TRY
    KRATOS_ERROR << "This function shouldn't be accessed here, use basic contact model instead."<<std::endl;
    KRATOS_CATCH("")
}

void DEM_parallel_bond::CalculateViscoDamping(double LocalRelVel[3],
                        double ViscoDampingLocalContactForce[3],
                        double indentation,
                        double equiv_visco_damp_coeff_normal,
                        double equiv_visco_damp_coeff_tangential,
                        bool& sliding,
                        int failure_id,
                        int i_neighbour_count,
                        SphericContinuumParticle* element1,
                        SphericContinuumParticle* element2) {
    KRATOS_TRY

    mUnbondedViscoDampingLocalContactForce[0] = 0.0;
    mUnbondedViscoDampingLocalContactForce[1] = 0.0;
    mUnbondedViscoDampingLocalContactForce[2] = 0.0;
    mBondedViscoDampingLocalContactForce[0] = 0.0;
    mBondedViscoDampingLocalContactForce[1] = 0.0;
    mBondedViscoDampingLocalContactForce[2] = 0.0;

    if (indentation > 0) {
        CalculateUnbondedViscoDampingForce(LocalRelVel, mUnbondedViscoDampingLocalContactForce, element1, element2);
    }

    if (!failure_id) {
        mBondedViscoDampingLocalContactForce[0] = -equiv_visco_damp_coeff_tangential * LocalRelVel[0];
        mBondedViscoDampingLocalContactForce[1] = -equiv_visco_damp_coeff_tangential * LocalRelVel[1];
        mBondedViscoDampingLocalContactForce[2] = -equiv_visco_damp_coeff_normal * LocalRelVel[2];
    }

    ViscoDampingLocalContactForce[0] = mUnbondedViscoDampingLocalContactForce[0] + mBondedViscoDampingLocalContactForce[0];
    ViscoDampingLocalContactForce[1] = mUnbondedViscoDampingLocalContactForce[1] + mBondedViscoDampingLocalContactForce[1];
    ViscoDampingLocalContactForce[2] = mUnbondedViscoDampingLocalContactForce[2] + mBondedViscoDampingLocalContactForce[2];

    double unbonded_normal_contact_force = mUnbondedLocalElasticContactForce2 + mUnbondedViscoDampingLocalContactForce[2];

    if (unbonded_normal_contact_force < 0.0) {
        mUnbondedViscoDampingLocalContactForce[2] = -1.0 * mUnbondedLocalElasticContactForce2;
        ViscoDampingLocalContactForce[2] = mUnbondedViscoDampingLocalContactForce[2] + mBondedViscoDampingLocalContactForce[2];
    }

    KRATOS_CATCH("") 
}

void DEM_parallel_bond::CalculateTangentialForces(double OldLocalElasticContactForce[3],
                        double LocalElasticContactForce[3],
                        double LocalElasticExtraContactForce[3],
                        double ViscoDampingLocalContactForce[3],
                        double LocalCoordSystem[3][3],
                        double LocalDeltDisp[3],
                        double LocalRelVel[3],
                        const double kt_el,
                        const double equiv_shear,
                        double& contact_tau,
                        double indentation,
                        double calculation_area,
                        double& failure_criterion_state,
                        SphericContinuumParticle* element1,
                        SphericContinuumParticle* element2,
                        int i_neighbour_count,
                        bool& sliding,
                        const ProcessInfo& r_process_info) {

    KRATOS_TRY

    // TODO: add [BOND_TAU_ZERO_DEVIATION]
    //const double& bond_tau_zero = (*mpProperties)[BOND_TAU_ZERO];
    //const double& internal_friction = (*mpProperties)[BOND_INTERNAL_FRICC];
    int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
    double current_tangential_force_module = 0.0;

    // The [mBondedScalingFactor] is divided into two direction [0] and [1]. June, 2022
    double BondedLocalElasticContactForce[2] = {0.0};

    // bond force
    if (!failure_type) {
        mAccumulatedBondedTangentialLocalDisplacement[0] += LocalDeltDisp[0];
        mAccumulatedBondedTangentialLocalDisplacement[1] += LocalDeltDisp[1];
        BondedLocalElasticContactForce[0] -= kt_el * mAccumulatedBondedTangentialLocalDisplacement[0]; // 0: first tangential
        BondedLocalElasticContactForce[1] -= kt_el * mAccumulatedBondedTangentialLocalDisplacement[1]; // 1: second tangential
    } else {
        //TODO: maybe a friction force due to the broekn bond should be added here
        BondedLocalElasticContactForce[0] = 0.0; // 0: first tangential
        BondedLocalElasticContactForce[1] = 0.0; // 1: second tangential
    }

    current_tangential_force_module = sqrt(BondedLocalElasticContactForce[0] * BondedLocalElasticContactForce[0]
                                                 + BondedLocalElasticContactForce[1] * BondedLocalElasticContactForce[1]);


    if (calculation_area){
        contact_tau = current_tangential_force_module / calculation_area;
    }

    //unbonded force
    double OldUnbondedLocalElasticContactForce[2] = {0.0};
    double UnbondedLocalElasticContactForce[2] = {0.0};
    double ActualTotalShearForce = 0.0;
    double max_admissible_shear_force = 0.0;
    double fraction = 0.0;

    if (indentation <= 0.0) {
        UnbondedLocalElasticContactForce[0] = 0.0;
        UnbondedLocalElasticContactForce[1] = 0.0;
    } else {
        // Here, unBondedScalingFactor[] = 1 - mBondedScalingFactor[]
        OldUnbondedLocalElasticContactForce[0] = (1 - mBondedScalingFactor[0]) * OldLocalElasticContactForce[0];
        OldUnbondedLocalElasticContactForce[1] = (1 - mBondedScalingFactor[1]) * OldLocalElasticContactForce[1];

        UnbondedLocalElasticContactForce[0] = OldUnbondedLocalElasticContactForce[0] - mKt * LocalDeltDisp[0];
        UnbondedLocalElasticContactForce[1] = OldUnbondedLocalElasticContactForce[1] - mKt * LocalDeltDisp[1];

        const double& equiv_tg_of_static_fri_ang = (*mpProperties)[STATIC_FRICTION];
        const double& equiv_tg_of_dynamic_fri_ang = (*mpProperties)[DYNAMIC_FRICTION];
        const double& equiv_friction_decay_coefficient = (*mpProperties)[FRICTION_DECAY];

        const double ShearRelVel = sqrt(LocalRelVel[0] * LocalRelVel[0] + LocalRelVel[1] * LocalRelVel[1]);
        // TODO: where this equation from? [equiv_friction]
        double equiv_friction = equiv_tg_of_dynamic_fri_ang + (equiv_tg_of_static_fri_ang - equiv_tg_of_dynamic_fri_ang) * exp(-equiv_friction_decay_coefficient * ShearRelVel);

        max_admissible_shear_force = (mUnbondedLocalElasticContactForce2 + mUnbondedViscoDampingLocalContactForce[2]) * equiv_friction;

        if (equiv_tg_of_static_fri_ang < 0.0 || equiv_tg_of_dynamic_fri_ang < 0.0) {
            KRATOS_ERROR << "The averaged friction is negative for one contact of element with Id: "<< element1->Id()<<std::endl;
        }

        const double tangential_contact_force_0 = UnbondedLocalElasticContactForce[0] + mUnbondedViscoDampingLocalContactForce[0];
        const double tangential_contact_force_1 = UnbondedLocalElasticContactForce[1] + mUnbondedViscoDampingLocalContactForce[1];

        ActualTotalShearForce = sqrt(tangential_contact_force_0 * tangential_contact_force_0 
                                    + tangential_contact_force_1 * tangential_contact_force_1);

        if (ActualTotalShearForce > max_admissible_shear_force) {
            const double ActualElasticShearForce = sqrt(UnbondedLocalElasticContactForce[0] * UnbondedLocalElasticContactForce[0] 
                                                        + UnbondedLocalElasticContactForce[1] * UnbondedLocalElasticContactForce[1]);

            const double dot_product = UnbondedLocalElasticContactForce[0] * mUnbondedViscoDampingLocalContactForce[0] 
                                        + UnbondedLocalElasticContactForce[1] * mUnbondedViscoDampingLocalContactForce[1];
            const double ViscoDampingLocalContactForceModule = sqrt(mUnbondedViscoDampingLocalContactForce[0] * mUnbondedViscoDampingLocalContactForce[0] +
                                                                    mUnbondedViscoDampingLocalContactForce[1] * mUnbondedViscoDampingLocalContactForce[1]);

            if (dot_product >= 0.0) {
                if (ActualElasticShearForce > max_admissible_shear_force) {
                    // if the ActualElasticShearForce is too big, it should be reduced
                    if(ActualElasticShearForce){
                        fraction = max_admissible_shear_force / ActualElasticShearForce;
                    }
                    UnbondedLocalElasticContactForce[0]      = UnbondedLocalElasticContactForce[0] * fraction;
                    UnbondedLocalElasticContactForce[1]      = UnbondedLocalElasticContactForce[1] * fraction;
                    mUnbondedViscoDampingLocalContactForce[0] = 0.0;
                    mUnbondedViscoDampingLocalContactForce[1] = 0.0;
                }
                else {
                    const double ActualViscousShearForce = max_admissible_shear_force - ActualElasticShearForce;
                    if(ViscoDampingLocalContactForceModule){
                       fraction = ActualViscousShearForce / ViscoDampingLocalContactForceModule;
                    }
                    mUnbondedViscoDampingLocalContactForce[0] *= fraction;
                    mUnbondedViscoDampingLocalContactForce[1] *= fraction;
                }
            }
            else {
                if (ViscoDampingLocalContactForceModule >= ActualElasticShearForce) {
                    if(ViscoDampingLocalContactForceModule){
                        fraction = (max_admissible_shear_force + ActualElasticShearForce) / ViscoDampingLocalContactForceModule;
                    }
                    mUnbondedViscoDampingLocalContactForce[0] *= fraction;
                    mUnbondedViscoDampingLocalContactForce[1] *= fraction;
                }
                else {
                    if(ActualElasticShearForce){
                        fraction = max_admissible_shear_force / ActualElasticShearForce;
                    }
                    UnbondedLocalElasticContactForce[0]      = UnbondedLocalElasticContactForce[0] * fraction;
                    UnbondedLocalElasticContactForce[1]      = UnbondedLocalElasticContactForce[1] * fraction;
                    mUnbondedViscoDampingLocalContactForce[0] = 0.0;
                    mUnbondedViscoDampingLocalContactForce[1] = 0.0;
                }
            }
            ViscoDampingLocalContactForce[0] = mBondedViscoDampingLocalContactForce[0] + mUnbondedViscoDampingLocalContactForce[0];
            ViscoDampingLocalContactForce[1] = mBondedViscoDampingLocalContactForce[1] + mUnbondedViscoDampingLocalContactForce[1];
            sliding = true;
        }
    }

    LocalElasticContactForce[0] = BondedLocalElasticContactForce[0] + UnbondedLocalElasticContactForce[0];
    LocalElasticContactForce[1] = BondedLocalElasticContactForce[1] + UnbondedLocalElasticContactForce[1];

    // Here, we only calculate the BondedScalingFactor and [unBondedScalingFactor = 1 - BondedScalingFactor].
    if (LocalElasticContactForce[0] && LocalElasticContactForce[1]) {
        mBondedScalingFactor[0] = BondedLocalElasticContactForce[0] / LocalElasticContactForce[0]; 
        mBondedScalingFactor[1] = BondedLocalElasticContactForce[1] / LocalElasticContactForce[1];
    } else {
        mBondedScalingFactor[0] = mBondedScalingFactor[1] = 0.0;
    }

    //for debug
    if (mDebugPrintingOption) {

        const long unsigned int& sphere_id = (*mpProperties)[DEBUG_PRINTING_ID_1];
        const long unsigned int& neigh_sphere_id = (*mpProperties)[DEBUG_PRINTING_ID_2];
        const double AccumulatedBondedTangentialLocalDisplacementModulus = sqrt(mAccumulatedBondedTangentialLocalDisplacement[0]*mAccumulatedBondedTangentialLocalDisplacement[0] + mAccumulatedBondedTangentialLocalDisplacement[1]*mAccumulatedBondedTangentialLocalDisplacement[1]);

        if ((element1->Id() == sphere_id) && (element2->Id() == neigh_sphere_id)) {
            std::ofstream normal_forces_file("delta_stress_tangential.txt", std::ios_base::out | std::ios_base::app);
            normal_forces_file << r_process_info[TIME] << " " << AccumulatedBondedTangentialLocalDisplacementModulus/*4*/ << " " << contact_tau/*5*/ << '\n'; 
            normal_forces_file.flush();
            normal_forces_file.close();
        }
    }

    KRATOS_CATCH("")
}

//*************************************
// Moment calculation
//*************************************

void DEM_parallel_bond::CalculateMoments(SphericContinuumParticle* element, 
                    SphericContinuumParticle* neighbor, 
                    double equiv_young, 
                    double distance, 
                    double calculation_area,
                    double LocalCoordSystem[3][3], 
                    double ElasticLocalRotationalMoment[3], 
                    double ViscoLocalRotationalMoment[3], 
                    double equiv_poisson, 
                    double indentation, 
                    double LocalElasticContactForce[3],
                    double normalLocalContactForce,
                    double GlobalElasticContactForces[3],
                    double LocalCoordSystem_2[3],
                    const int i_neighbor_count) 
{
    KRATOS_TRY

    int failure_type = element->mIniNeighbourFailureId[i_neighbor_count];

    if (failure_type == 0) {
            ComputeParticleRotationalMoments(element, 
                                    neighbor, 
                                    equiv_young, 
                                    distance, 
                                    calculation_area,
                                    LocalCoordSystem, 
                                    ElasticLocalRotationalMoment, 
                                    ViscoLocalRotationalMoment, 
                                    equiv_poisson, 
                                    indentation, 
                                    LocalElasticContactForce);
    }
                            
    double LocalUnbondElasticContactForce[3] = {0.0};
    double GlobalUnbondElasticContactForce[3] = {0.0};
    LocalUnbondElasticContactForce[0] = (1 - mBondedScalingFactor[0]) * LocalElasticContactForce[0];
    LocalUnbondElasticContactForce[1] = (1 - mBondedScalingFactor[1]) * LocalElasticContactForce[1];
    LocalUnbondElasticContactForce[2] = (1 - mBondedScalingFactor[2]) * LocalElasticContactForce[2];
    GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalUnbondElasticContactForce, GlobalUnbondElasticContactForce);

    DemContact::ComputeParticleContactMoments(normalLocalContactForce,
                                             GlobalUnbondElasticContactForce,
                                             LocalCoordSystem_2,
                                             element,
                                             neighbor,
                                             indentation,
                                             i_neighbor_count);

    KRATOS_CATCH("")
}

double DEM_parallel_bond::GetYoungModulusForComputingRotationalMoments(const double& equiv_young){
        const double bond_equiv_young = (*mpProperties)[BOND_YOUNG_MODULUS];
        return bond_equiv_young;
    }

void DEM_parallel_bond::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                SphericContinuumParticle* neighbor,
                                                double equiv_young,
                                                double distance,
                                                double calculation_area,
                                                double LocalCoordSystem[3][3],
                                                double ElasticLocalRotationalMoment[3],
                                                double ViscoLocalRotationalMoment[3],
                                                double equiv_poisson,
                                                double indentation,
                                                double LocalElasticContactForce[3]) {

    KRATOS_TRY

    double LocalDeltaRotatedAngle[3]    = {0.0};
    double LocalDeltaAngularVelocity[3] = {0.0};

    array_1d<double, 3> GlobalDeltaRotatedAngle;
    noalias(GlobalDeltaRotatedAngle) = element->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE) 
                                        - neighbor->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
    array_1d<double, 3> GlobalDeltaAngularVelocity;
    noalias(GlobalDeltaAngularVelocity) = element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY) 
                                        - neighbor->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

    GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaRotatedAngle, LocalDeltaRotatedAngle);
    GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalDeltaAngularVelocity, LocalDeltaAngularVelocity);

    const double equivalent_radius = std::sqrt(calculation_area / Globals::Pi);
    const double element_mass  = element->GetMass();
    const double neighbor_mass = neighbor->GetMass();
    const double equiv_mass    = element_mass * neighbor_mass / (element_mass + neighbor_mass);

    const double bond_equiv_young = GetYoungModulusForComputingRotationalMoments(equiv_young);
    
    double kn_el = bond_equiv_young * calculation_area / distance;
    double kt_el = kn_el / (*mpProperties)[BOND_KNKS_RATIO];

    const double Inertia_I     = 0.25 * Globals::Pi * equivalent_radius * equivalent_radius * equivalent_radius * equivalent_radius;
    const double Inertia_J     = 2.0 * Inertia_I; // This is the polar inertia

    const double& damping_gamma = (*mpProperties)[DAMPING_GAMMA];

    //Viscous parameter taken from Olmedo et al., 'Discrete element model of the dynamic response of fresh wood stems to impact'
    array_1d<double, 3> visc_param;
    visc_param[0] = 2.0 * damping_gamma * std::sqrt(equiv_mass * bond_equiv_young * Inertia_I / distance); // OLMEDO
    visc_param[1] = 2.0 * damping_gamma * std::sqrt(equiv_mass * bond_equiv_young * Inertia_I / distance); // OLMEDO
    visc_param[2] = 2.0 * damping_gamma * std::sqrt(equiv_mass * bond_equiv_young * Inertia_J / distance); // OLMEDO

    double aux = 0.0;
    aux = (element->GetRadius() + neighbor->GetRadius()) / distance; // This is necessary because if spheres are not tangent the DeltaAngularVelocity has to be interpolated
   
 
    array_1d<double, 3> LocalEffDeltaRotatedAngle;
    LocalEffDeltaRotatedAngle[0] = LocalDeltaRotatedAngle[0] * aux;
    LocalEffDeltaRotatedAngle[1] = LocalDeltaRotatedAngle[1] * aux;
    LocalEffDeltaRotatedAngle[2] = LocalDeltaRotatedAngle[2] * aux;

    array_1d<double, 3> LocalEffDeltaAngularVelocity;
    LocalEffDeltaAngularVelocity[0] = LocalDeltaAngularVelocity[0] * aux;
    LocalEffDeltaAngularVelocity[1] = LocalDeltaAngularVelocity[1] * aux;
    LocalEffDeltaAngularVelocity[2] = LocalDeltaAngularVelocity[2] * aux;

    ElasticLocalRotationalMoment[0] = -kn_el / calculation_area * Inertia_I * LocalEffDeltaRotatedAngle[0];
    ElasticLocalRotationalMoment[1] = -kn_el / calculation_area * Inertia_I * LocalEffDeltaRotatedAngle[1];
    ElasticLocalRotationalMoment[2] = -kt_el / calculation_area * Inertia_J * LocalEffDeltaRotatedAngle[2];

    // Bond rotational 'friction' based on particle rolling fricton 
    //Not damping but simple implementation to help energy dissipation
    
    double LocalElement1AngularVelocity[3] = {0.0};
    array_1d<double, 3> GlobalElement1AngularVelocity;
    noalias(GlobalElement1AngularVelocity) = element->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
    GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalElement1AngularVelocity, LocalElement1AngularVelocity);
    double element1AngularVelocity_modulus = sqrt(LocalElement1AngularVelocity[0] * LocalElement1AngularVelocity[0] + 
                                                LocalElement1AngularVelocity[1] * LocalElement1AngularVelocity[1] +
                                                LocalElement1AngularVelocity[2] * LocalElement1AngularVelocity[2]);

    if (element1AngularVelocity_modulus){
        array_1d<double, 3> other_to_me_vect;
        noalias(other_to_me_vect) = element->GetGeometry()[0].Coordinates() - neighbor->GetGeometry()[0].Coordinates();
        double bond_center_point_to_element1_mass_center_distance = DEM_MODULUS_3(other_to_me_vect) / 2; //Here, this only works for sphere particles
    
        array_1d<double, 3> element1AngularVelocity_normalise;
        element1AngularVelocity_normalise[0] = LocalElement1AngularVelocity[0] / element1AngularVelocity_modulus;
        element1AngularVelocity_normalise[1] = LocalElement1AngularVelocity[1] / element1AngularVelocity_modulus;
        element1AngularVelocity_normalise[2] = LocalElement1AngularVelocity[2] / element1AngularVelocity_modulus;

        Properties& properties_of_this_contact = element->GetProperties().GetSubProperties(neighbor->GetProperties().Id());
        
        double BondedLocalElasticContactForce[3] = {0.0};
        BondedLocalElasticContactForce[0] = mBondedScalingFactor[0] * LocalElasticContactForce[0];
        BondedLocalElasticContactForce[1] = mBondedScalingFactor[1] * LocalElasticContactForce[1];
        BondedLocalElasticContactForce[2] = mBondedScalingFactor[2] * LocalElasticContactForce[2];

        ViscoLocalRotationalMoment[0] = - element1AngularVelocity_normalise[0] * std::abs(BondedLocalElasticContactForce[2]) * bond_center_point_to_element1_mass_center_distance 
                                        * properties_of_this_contact[ROLLING_FRICTION]; 

        ViscoLocalRotationalMoment[1] = - element1AngularVelocity_normalise[1] * std::abs(BondedLocalElasticContactForce[2]) * bond_center_point_to_element1_mass_center_distance 
                                        * properties_of_this_contact[ROLLING_FRICTION]; 

        ViscoLocalRotationalMoment[2] = - element1AngularVelocity_normalise[2] * std::abs(BondedLocalElasticContactForce[2]) * bond_center_point_to_element1_mass_center_distance 
                                        * properties_of_this_contact[ROLLING_FRICTION]; 

    } else {
        ViscoLocalRotationalMoment[0] = 0.0;
        ViscoLocalRotationalMoment[1] = 0.0;
        ViscoLocalRotationalMoment[2] = 0.0;
    }
     
    KRATOS_CATCH("")
}//ComputeParticleRotationalMoments


//*************************************
// Bond failure checking
//*************************************

void DEM_parallel_bond::CheckFailure(const int i_neighbour_count, 
                                        SphericContinuumParticle* element1, 
                                        SphericContinuumParticle* element2,
                                        double& contact_sigma,
                                        double& contact_tau, 
                                        double LocalElasticContactForce[3],
                                        double ViscoDampingLocalContactForce[3],
                                        double ElasticLocalRotationalMoment[3],
                                        double ViscoLocalRotationalMoment[3]){
    
    KRATOS_TRY

    int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

    if (failure_type == 0) {

        //parameters
        const double& bond_sigma_max = (*mpProperties)[BOND_SIGMA_MAX];
        //const double& bond_sigma_max_deviation = (*mpProperties)[BOND_SIGMA_MAX_DEVIATION];
        const double& bond_tau_zero = (*mpProperties)[BOND_TAU_ZERO];
        //const double& bond_tau_zero_deviation = (*mpProperties)[BOND_TAU_ZERO_DEVIATION];
        const double& bond_interanl_friction = (*mpProperties)[BOND_INTERNAL_FRICC];
        const double& bond_rotational_moment_coefficient_normal =(*mpProperties)[BOND_ROTATIONAL_MOMENT_COEFFICIENT_NORMAL];
        const double& bond_rotational_moment_coefficient_tangential =(*mpProperties)[BOND_ROTATIONAL_MOMENT_COEFFICIENT_TANGENTIAL];
        
        double bond_rotational_moment[3] = {0.0};
        bond_rotational_moment[0]  = ElasticLocalRotationalMoment[0];
        bond_rotational_moment[1]  = ElasticLocalRotationalMoment[1];
        bond_rotational_moment[2]  = ElasticLocalRotationalMoment[2];
        double bond_rotational_moment_normal_modulus = 0.0;
        double bond_rotational_moment_tangential_modulus = 0.0;

        bond_rotational_moment_normal_modulus     = std::abs(bond_rotational_moment[2]);
        bond_rotational_moment_tangential_modulus = sqrt(bond_rotational_moment[0] * bond_rotational_moment[0]
                                                    + bond_rotational_moment[1] * bond_rotational_moment[1]);

        const double my_radius         = element1->GetRadius();
        const double other_radius      = element2->GetRadius();
        const double bond_radius_factor = (*mpProperties)[BOND_RADIUS_FACTOR];
        double bond_radius = std::min(my_radius, other_radius) * bond_radius_factor;

        const double I = 0.25 * Globals::Pi * bond_radius * bond_radius * bond_radius * bond_radius;
        const double J = 2.0 * I; // This is the polar inertia

        double bond_current_tau_max = bond_tau_zero;

        if (contact_sigma >= 0) {
            bond_current_tau_max += tan(bond_interanl_friction * Globals::Pi / 180.0) * contact_sigma;
        }

        if (contact_sigma < 0.0  /*break only in tension*/
                && (-1 * contact_sigma + bond_rotational_moment_coefficient_normal * bond_rotational_moment_tangential_modulus * bond_radius / I > bond_sigma_max) 
                && !(*mpProperties)[IS_UNBREAKABLE]) 
        { //for normal
            failure_type = 4; // failure in tension
            contact_sigma = 0.0;
            contact_tau = 0.0;
            LocalElasticContactForce[0] *= (1 - mBondedScalingFactor[0]);      
            LocalElasticContactForce[1] *= (1 - mBondedScalingFactor[1]);      
            LocalElasticContactForce[2]  = mUnbondedLocalElasticContactForce2;
            ViscoDampingLocalContactForce[0] = mUnbondedViscoDampingLocalContactForce[0];
            ViscoDampingLocalContactForce[1] = mUnbondedViscoDampingLocalContactForce[1];
            ViscoDampingLocalContactForce[2] = mUnbondedViscoDampingLocalContactForce[2];
            ElasticLocalRotationalMoment[0] = 0.0;
            ElasticLocalRotationalMoment[1] = 0.0;
            ElasticLocalRotationalMoment[2] = 0.0;
            ViscoLocalRotationalMoment[0] = 0.0;
            ViscoLocalRotationalMoment[1] = 0.0;
            ViscoLocalRotationalMoment[2] = 0.0;
        }
        else if(( std::abs(contact_tau) + bond_rotational_moment_coefficient_tangential * bond_rotational_moment_normal_modulus * bond_radius / J > bond_current_tau_max) 
            && !(*mpProperties)[IS_UNBREAKABLE]) 
        { //for tangential 
            failure_type = 2; // failure in shear
            contact_sigma = 0.0;
            contact_tau = 0.0;
            //If bond break in shear, the normal compressive force should still be there like before
            LocalElasticContactForce[0] *= (1 - mBondedScalingFactor[0]);      
            LocalElasticContactForce[1] *= (1 - mBondedScalingFactor[1]);      
            //LocalElasticContactForce[2]  = mUnbondedLocalElasticContactForce2; 
            ViscoDampingLocalContactForce[0] = mUnbondedViscoDampingLocalContactForce[0];
            ViscoDampingLocalContactForce[1] = mUnbondedViscoDampingLocalContactForce[1];
            //ViscoDampingLocalContactForce[2] = mUnbondedViscoDampingLocalContactForce[2];
            ElasticLocalRotationalMoment[0] = 0.0;
            ElasticLocalRotationalMoment[1] = 0.0;
            ElasticLocalRotationalMoment[2] = 0.0;
            ViscoLocalRotationalMoment[0] = 0.0;
            ViscoLocalRotationalMoment[1] = 0.0;
            ViscoLocalRotationalMoment[2] = 0.0;
        } 
    }

    KRATOS_CATCH("")    
}//CheckFailure
    
} //namespace Kratos