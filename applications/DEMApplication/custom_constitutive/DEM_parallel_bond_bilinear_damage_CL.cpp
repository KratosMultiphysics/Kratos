/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu
// Date: Oct 2023
/////////////////////////////////////////////////

// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <algorithm>

// Project includes
#include "DEM_parallel_bond_bilinear_damage_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos{

DEMContinuumConstitutiveLaw::Pointer DEM_parallel_bond_bilinear_damage::Clone() const{
    DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_parallel_bond_bilinear_damage(*this));
    return p_clone;
}

//*************************************
// Parameters preparation
//*************************************

void DEM_parallel_bond_bilinear_damage::Initialize(SphericContinuumParticle* element1,
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

void DEM_parallel_bond_bilinear_damage::Check(Properties::Pointer pProp) const {

    BaseClassType::Check(pProp);

    if (!pProp->Has(FRACTURE_ENERGY_NORMAL)) {
        KRATOS_WARNING("DEM")<<"\nWARNING: Variable FRACTURE_ENERGY_NORMAL was not found in cthe Properties when using DEM_parallel_bond_bilinear_damage. A default value of 0.0 was assigned.\n\n";
        pProp->GetValue(FRACTURE_ENERGY_NORMAL) = 0.0;
    }

    if (!pProp->Has(FRACTURE_ENERGY_TANGENTIAL)) {
        KRATOS_WARNING("DEM")<<"\nWARNING: Variable FRACTURE_ENERGY_TANGENTIAL was not found in cthe Properties when using DEM_parallel_bond_bilinear_damage. A default value of 0.0 was assigned.\n\n";
        pProp->GetValue(FRACTURE_ENERGY_TANGENTIAL) = 0.0;
    }
}

//*************************************
// Force calculation
//*************************************

double DEM_parallel_bond_bilinear_damage::ComputeNormalUnbondedForce(double indentation){
    
    KRATOS_TRY

    KRATOS_ERROR << "This function shouldn't be accessed here, use basic contact model instead."<<std::endl
                << "Maybe you are using \"DEM_parallel_bond_bilinear_damage\" for DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME." <<std::endl
                << "Unfortunately, you can only input one of the names listed below." <<std::endl
                << "1. DEM_parallel_bond_bilinear_damage_Linear" <<std::endl 
                << "2. DEM_parallel_bond_bilinear_damage_Hertz" <<std::endl
                << "3. DEM_parallel_bond_bilinear_damage_Quadratic" <<std::endl;
    
    KRATOS_CATCH("")
}

void DEM_parallel_bond_bilinear_damage::CalculateForces(const ProcessInfo& r_process_info,
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
                            double indentation_particle,
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

    //
    //In this CL, we calculate the normal and tangential forces together in this function for calculting a uniform damage in normal and tangential direction.
    //
    
    //*******************Calculate NormalForces (and also the tangential bond force)************************
    int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
    const double bonded_indentation = indentation - mInitialIndentationForBondedPart; 
                                                                                    
    double BondedLocalElasticContactForce2 = 0.0; //normal forces
    // The [mBondedScalingFactor] is divided into two direction [0] and [1]. June, 2022
    double OldBondedLocalElasticContactForce[2] = {0.0};
    OldBondedLocalElasticContactForce[0] = mBondedScalingFactor[0] * OldLocalElasticContactForce[0];
    OldBondedLocalElasticContactForce[1] = mBondedScalingFactor[1] * OldLocalElasticContactForce[1];
    double BondedLocalElasticContactForce[2] = {0.0}; 
    double current_tangential_force_module = 0.0;

    const double bond_sigma_max = (*mpProperties)[BOND_SIGMA_MAX]; //tension limit
    const double bond_tau_zero = (*mpProperties)[BOND_TAU_ZERO];
    const double internal_friction = (*mpProperties)[BOND_INTERNAL_FRICC];
    const double fracture_energy_normal = (*mpProperties)[FRACTURE_ENERGY_NORMAL];
    const double fracture_energy_tangential = (*mpProperties)[FRACTURE_ENERGY_TANGENTIAL];
    const double delta_at_undamaged_peak_normal = bond_sigma_max * calculation_area / kn_el;
    const double delta_at_failure_point_normal = (2.0 * fracture_energy_normal) / bond_sigma_max;

    /*if (std::abs(delta_at_failure_point_normal) < std::abs(delta_at_undamaged_peak_normal) && fracture_energy_normal){
        double suggested_fracture_energy = delta_at_undamaged_peak_normal * bond_sigma_max / 2.0;
        KRATOS_INFO("DEM") << "The [Fracture_energy_normal] is too small! It shouble be bigger than " << suggested_fracture_energy << std::endl;
        KRATOS_ERROR<< "The [Fracture_energy_normal] is too small!" << std::endl;
    }*/
    
  //**********************************************
    //normal
    const double initial_limit_force = bond_sigma_max * calculation_area;
    double DamageEnergyCoeffNormal = 0.0;

    if (bond_sigma_max) {
        DamageEnergyCoeffNormal = 2.0 * fracture_energy_normal * kn_el / (calculation_area * bond_sigma_max * bond_sigma_max) - 1.0;
    } else {
        DamageEnergyCoeffNormal = 0.0;
    }

    KRATOS_ERROR_IF(DamageEnergyCoeffNormal > 30.0) << "Damage energy is too big!" << std::endl;

    if (DamageEnergyCoeffNormal < 0.0) {
        DamageEnergyCoeffNormal = 0.0;
    }

    double k_softening = 0.0;
    double limit_force = 0.0;

    if (DamageEnergyCoeffNormal) {
        k_softening = kn_el / DamageEnergyCoeffNormal;
    }

    double kn_updated = (1.0 - mDamageNormal) * kn_el;
    
    double delta_accumulated = 0.0;

    double returned_by_mapping_force = 0.0;

    //double BondedLocalElasticContactForce2 = 0.0;
    //const double bonded_indentation = indentation - mInitialIndentationForBondedPart;
    double current_normal_force_module = 0.0;

    if (bonded_indentation >= 0.0) { //COMPRESSION
        if (!failure_type) {
            BondedLocalElasticContactForce2 = kn_updated * bonded_indentation;
            delta_accumulated = bonded_indentation;
        } else {
            BondedLocalElasticContactForce2 = 0.0;
        }
    } else { //tension

        if (!failure_type) {
            if (DamageEnergyCoeffNormal) {
                limit_force = initial_limit_force * (1.0 + k_softening / kn_el) * kn_updated / (kn_updated + k_softening);
            } else {
                limit_force = initial_limit_force;
            }

            current_normal_force_module = fabs(kn_updated * bonded_indentation);
            
            BondedLocalElasticContactForce2 = kn_updated * bonded_indentation;

            delta_accumulated = current_normal_force_module / kn_updated;

            returned_by_mapping_force = current_normal_force_module;

            if ((current_normal_force_module > limit_force) && !(*mpProperties)[IS_UNBREAKABLE]) {

                if (DamageEnergyCoeffNormal) { // the material can sustain further damage, not failure yet

                    const double delta_at_undamaged_peak = initial_limit_force / kn_el;

                    returned_by_mapping_force = initial_limit_force - k_softening * (delta_accumulated - delta_at_undamaged_peak);

                    if (returned_by_mapping_force < 0.0) {
                        returned_by_mapping_force = 0.0;
                    }

                    BondedLocalElasticContactForce2 = -returned_by_mapping_force;

                    mDamageNormal = 1.0 - (returned_by_mapping_force / delta_accumulated) / kn_el;

                    if (mDamageNormal > mDamageThresholdTolerance) {
                        failure_type = 4; // failure by traction
                        BondedLocalElasticContactForce2 = 0.0;
                        mDamageNormal = 1.0;
                    }
                } else { // Fully fragile behaviour

                    failure_type = 4; // failure by traction
                    BondedLocalElasticContactForce2 = 0.0;
                    mDamageNormal = 1.0;
                }
            }
        } else {
            BondedLocalElasticContactForce2 = 0.0;
        }
    }

    //tangential bond force
    k_softening = 0.0;
    double tau_strength = 0.0;

    //double OldBondedLocalElasticContactForce[2] = {0.0};
    //OldBondedLocalElasticContactForce[0] = mBondedScalingFactor * OldLocalElasticContactForce[0];
    //OldBondedLocalElasticContactForce[1] = mBondedScalingFactor * OldLocalElasticContactForce[1];
    double DamageEnergyCoeffTangential = 0.0;

    if (bond_sigma_max) {
        DamageEnergyCoeffTangential = 2.0 * fracture_energy_tangential * kt_el / (calculation_area * bond_tau_zero * bond_tau_zero) - 1.0;
    } else {
        DamageEnergyCoeffTangential = 0.0;
    }

    KRATOS_ERROR_IF(DamageEnergyCoeffTangential > 30.0) << "Damage energy is too big!" << std::endl;

    if (DamageEnergyCoeffTangential < 0.0) {
        DamageEnergyCoeffTangential = 0.0;
    }

    if (DamageEnergyCoeffTangential) {
        k_softening = kt_el / DamageEnergyCoeffTangential;
    }

    //int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
    double kt_updated = (1.0 - mDamageTangential) * kt_el;

    //double BondedLocalElasticContactForce[2] = {0.0};

    if (!failure_type) {
        BondedLocalElasticContactForce[0] = OldBondedLocalElasticContactForce[0] - kt_updated * LocalDeltDisp[0]; // 0: first tangential
        BondedLocalElasticContactForce[1] = OldBondedLocalElasticContactForce[1] - kt_updated * LocalDeltDisp[1]; // 1: second tangential
    } else {
        BondedLocalElasticContactForce[0] = 0.0; // 0: first tangential
        BondedLocalElasticContactForce[1] = 0.0; // 1: second tangential
    }

    delta_accumulated = 0.0;
    //double current_tangential_force_module = 0.0;
    returned_by_mapping_force = 0.0;

    if (!failure_type) { // This means it has not broken yet

        current_tangential_force_module = sqrt(BondedLocalElasticContactForce[0] * BondedLocalElasticContactForce[0]
                                                + BondedLocalElasticContactForce[1] * BondedLocalElasticContactForce[1]);

        delta_accumulated = current_tangential_force_module / kt_updated;

        returned_by_mapping_force = current_tangential_force_module;

        contact_sigma = LocalElasticContactForce[2] / calculation_area;
        contact_tau = current_tangential_force_module / calculation_area;

        double updated_max_tau_strength = bond_tau_zero;

        if (contact_sigma >= 0) {
            updated_max_tau_strength += internal_friction * contact_sigma;
        }
        tau_strength = updated_max_tau_strength * (1.0 + k_softening / kt_el) * kt_updated / (kt_updated + k_softening);

        if ((contact_tau > tau_strength) && !(*mpProperties)[IS_UNBREAKABLE]) { // damage

            if (DamageEnergyCoeffTangential) { // the material can sustain further damage, not failure yet

                const double delta_at_undamaged_peak = updated_max_tau_strength * calculation_area / kt_el;

                delta_accumulated = current_tangential_force_module / kt_updated;

                returned_by_mapping_force = updated_max_tau_strength * calculation_area - k_softening * (delta_accumulated - delta_at_undamaged_peak);

                if (returned_by_mapping_force < 0.0) {
                    returned_by_mapping_force = 0.0;
                }

                if (current_tangential_force_module) {
                    BondedLocalElasticContactForce[0] = (returned_by_mapping_force / current_tangential_force_module) * BondedLocalElasticContactForce[0];
                    BondedLocalElasticContactForce[1] = (returned_by_mapping_force / current_tangential_force_module) * BondedLocalElasticContactForce[1];
                }

                mDamageTangential = 1.0 - (returned_by_mapping_force / delta_accumulated) / kt_el; // This is 1 - quotient of K_damaged/K_intact

                if (mDamageTangential > mDamageThresholdTolerance) {
                    failure_type = 2; // failure by shear
                    BondedLocalElasticContactForce[0] = 0.0;
                    BondedLocalElasticContactForce[1] = 0.0;
                    mDamageTangential = 1.0;
                }
            } else { // Fully fragile behaviour

                failure_type = 2; // failure by shear
                BondedLocalElasticContactForce[0] = 0.0;
                BondedLocalElasticContactForce[1] = 0.0;
                mDamageTangential = 1.0;
            }
        }
    }

    if (indentation_particle > 0.0) {
        mUnbondedLocalElasticContactForce2 = ComputeNormalUnbondedForce(indentation_particle);
    } else {
        mUnbondedLocalElasticContactForce2 = 0.0;
    }
    
    current_tangential_force_module = sqrt(BondedLocalElasticContactForce[0] * BondedLocalElasticContactForce[0]
                                                 + BondedLocalElasticContactForce[1] * BondedLocalElasticContactForce[1]);
        
    if(calculation_area){
        contact_sigma = BondedLocalElasticContactForce2 / calculation_area;
        contact_tau = current_tangential_force_module / calculation_area;
    }

    LocalElasticContactForce[2] = BondedLocalElasticContactForce2 + mUnbondedLocalElasticContactForce2;

    if (LocalElasticContactForce[2]) {
        mBondedScalingFactor[2] = BondedLocalElasticContactForce2 / LocalElasticContactForce[2]; 
    } else {
        mBondedScalingFactor[2] = 0.0;
    }       
    
    //********************Calculate ViscoDampingCoeff***************************
    if (mDamageReal < 1.0){
        const double my_mass    = element1->GetMass();
        const double other_mass = element2->GetMass();

        const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);

        // TODO: check the value of gamma
        const double equiv_gamma = (*mpProperties)[DAMPING_GAMMA];

        double kn_el_update = (1 - mDamageNormal) * kn_el;
        double kt_el_update = (1 - mDamageTangential) * kt_el;

        equiv_visco_damp_coeff_normal     = 2.0 * equiv_gamma * sqrt(equiv_mass * kn_el_update);
        equiv_visco_damp_coeff_tangential = 2.0 * equiv_gamma * sqrt(equiv_mass * kt_el_update);

    } else {
        equiv_visco_damp_coeff_normal     = 0.0;
        equiv_visco_damp_coeff_tangential = 0.0;
    }


    //************************Calculate ViscoDamping*****************************

    BaseClassType::CalculateViscoDamping(LocalRelVel,
                        ViscoDampingLocalContactForce,
                        indentation_particle,
                        equiv_visco_damp_coeff_normal,
                        equiv_visco_damp_coeff_tangential,
                        sliding,
                        element1->mIniNeighbourFailureId[i_neighbour_count],
                        i_neighbour_count,
                        element1,
                        element2);

    // Tangential forces are calculated after the viscodamping because the frictional limit bounds the sum of elastic plus viscous forces
    
    //************************Calculate tangential forces*****************************

    //unbonded force
    double OldUnbondedLocalElasticContactForce[2] = {0.0};
    double UnbondedLocalElasticContactForce[2] = {0.0};
    double ActualTotalShearForce = 0.0;
    double max_admissible_shear_force = 0.0;
    double fraction = 0.0;

    if (indentation_particle <= 0.0) {
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
    if (LocalElasticContactForce[0]) {
        mBondedScalingFactor[0] = BondedLocalElasticContactForce[0] / LocalElasticContactForce[0]; 
    } else {
        mBondedScalingFactor[0] = 0.0;
    }

    if (LocalElasticContactForce[1]) { 
        mBondedScalingFactor[1] = BondedLocalElasticContactForce[1] / LocalElasticContactForce[1];
    } else {
        mBondedScalingFactor[1] = 0.0;
    }

    //for debug
    if (mDebugPrintingOption) {

        const long unsigned int& sphere_id = (*mpProperties)[DEBUG_PRINTING_ID_1];
        const long unsigned int& neigh_sphere_id = (*mpProperties)[DEBUG_PRINTING_ID_2];
        const double AccumulatedBondedTangentialLocalDisplacementModulus = sqrt(mAccumulatedBondedTangentialLocalDisplacement[0]*mAccumulatedBondedTangentialLocalDisplacement[0] + mAccumulatedBondedTangentialLocalDisplacement[1]*mAccumulatedBondedTangentialLocalDisplacement[1]);

        if ((element1->Id() == sphere_id) && (element2->Id() == neigh_sphere_id)) {
            std::ofstream normal_forces_file("delta_stress.txt", std::ios_base::out | std::ios_base::app);
            normal_forces_file << r_process_info[TIME] << " " << bonded_indentation/*2*/ << " " << contact_sigma/*3*/ << " "
                                << AccumulatedBondedTangentialLocalDisplacementModulus/*4*/ << " " << contact_tau/*5*/ << " " 
                                << mDamageReal/*6*/ << " " << mDamageNormal << " " << bonded_indentation << " " << '\n'; 
            normal_forces_file.flush();
            normal_forces_file.close();
        }
    }

    if (mDamageNormal > mDamageTangential){
        mDamageTangential = mDamageNormal;
    } else {
        mDamageNormal = mDamageTangential;
    }

    KRATOS_CATCH("") 
}


//*************************************
// Moment calculation
//*************************************

void DEM_parallel_bond_bilinear_damage::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                SphericContinuumParticle* neighbor,
                                                double equiv_young,
                                                double distance,
                                                double calculation_area,
                                                double LocalCoordSystem[3][3],
                                                double ElasticLocalRotationalMoment[3],
                                                double ViscoLocalRotationalMoment[3],
                                                double equiv_poisson,
                                                double indentation,
                                                double LocalContactForce[3]) {

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
    
    double kn_el = (1 - mDamageNormal) * bond_equiv_young * calculation_area / distance;
    double kt_el = (1 - mDamageTangential) * kn_el / (*mpProperties)[BOND_KNKS_RATIO];

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
    /*
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
    */
     
    KRATOS_CATCH("")
}//ComputeParticleRotationalMoments


//*************************************
// Bond failure checking
//*************************************

void DEM_parallel_bond_bilinear_damage::CheckFailure(const int i_neighbour_count, 
                                        SphericContinuumParticle* element1, 
                                        SphericContinuumParticle* element2,
                                        double& contact_sigma,
                                        double& contact_tau, 
                                        double LocalElasticContactForce[3],
                                        double ViscoDampingLocalContactForce[3],
                                        double ElasticLocalRotationalMoment[3],
                                        double ViscoLocalRotationalMoment[3]){
    
    KRATOS_TRY

    //The momonets do not contribute to failure in current version

    KRATOS_CATCH("")    
}//CheckFailure

    
} //namespace Kratos