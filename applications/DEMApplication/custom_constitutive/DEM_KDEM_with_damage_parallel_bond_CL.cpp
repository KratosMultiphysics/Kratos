// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <cstdio>

// Project includes
#include "DEM_KDEM_with_damage_parallel_bond_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_with_damage_parallel_bond::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_with_damage_parallel_bond(*this));
        return p_clone;
    }

    void DEM_KDEM_with_damage_parallel_bond::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_with_damage_parallel_bond to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEM_KDEM_with_damage_parallel_bond::Check(Properties::Pointer pProp) const {

        DEM_KDEM_with_damage::Check(pProp);

        if (!pProp->Has(BONDED_MATERIAL_YOUNG_MODULUS)) {
            KRATOS_WARNING("DEM")<<"\nWARNING: Variable BONDED_MATERIAL_YOUNG_MODULUS was not found in the Properties when using DEM_KDEM_with_damage_parallel_bond. A default value of 0.0 was assigned.\n\n";
            pProp->GetValue(BONDED_MATERIAL_YOUNG_MODULUS) = 0.0;
        }
        if (!pProp->Has(FRACTURE_ENERGY)) {
            KRATOS_WARNING("DEM")<<"\nWARNING: Variable FRACTURE_ENERGY was not found in the Properties when using DEM_KDEM_with_damage_parallel_bond. A default value of 0.0 was assigned.\n\n";
            pProp->GetValue(FRACTURE_ENERGY) = 0.0;
        }
    }

    void DEM_KDEM_with_damage_parallel_bond::Initialize(SphericContinuumParticle* element) {

        mDebugPrintingOption = false;
        Properties::Pointer pProp = element->pGetProperties();
        if (!pProp->Has(DEBUG_PRINTING_OPTION)) {
            mDebugPrintingOption = false;
        } else {
            mDebugPrintingOption = bool(pProp->GetValue(DEBUG_PRINTING_OPTION));
        }
        if (mDebugPrintingOption) {
            if (!pProp->Has(DEBUG_PRINTING_ID_1) || !pProp->Has(DEBUG_PRINTING_ID_2)) {
                KRATOS_WARNING("DEM") << "\nWARNING: We are currently in DEBUG PRINTING mode, so the ids of the two particles involved must be given.\n\n";
            }
        }
    }

    void DEM_KDEM_with_damage_parallel_bond::CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                             double equiv_poisson, double calculation_area, SphericContinuumParticle* element1, SphericContinuumParticle* element2, double indentation) {

        KRATOS_TRY

        //TODO: Sometimes we do not compute mean values this way. Sometimes we use 2xy/(x+y)
        const double unbonded_equivalent_young = 0.5 * (element1->GetProperties()[YOUNG_MODULUS] + element2->GetProperties()[YOUNG_MODULUS]);
        const double unbonded_equivalent_shear = unbonded_equivalent_young / (2.0 * (1 + equiv_poisson));

        mUnbondedNormalElasticConstant = unbonded_equivalent_young * calculation_area / initial_dist;
        mUnbondedTangentialElasticConstant = unbonded_equivalent_shear * calculation_area / initial_dist;

        const double bonded_equiv_young = 0.5 * (element1->GetProperties()[BONDED_MATERIAL_YOUNG_MODULUS] + element2->GetProperties()[BONDED_MATERIAL_YOUNG_MODULUS]);
        const double bonded_equiv_shear = bonded_equiv_young / (2.0 * (1 + equiv_poisson));

        kn_el = bonded_equiv_young * calculation_area / initial_dist;
        kt_el = bonded_equiv_shear * calculation_area / initial_dist;

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage_parallel_bond::CalculateForces(const ProcessInfo& r_process_info,
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
                                double& accumulated_damage,
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
                accumulated_damage,
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

        FindMaximumValueOfNormalAndTangentialDamageComponents(element1, element2);

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage_parallel_bond::FindMaximumValueOfNormalAndTangentialDamageComponents(SphericContinuumParticle* element1,
                                                                                                   SphericContinuumParticle* element2) {

        KRATOS_TRY

        mDamageNormal = std::max(mDamageNormal, mDamageTangential);
        mDamageTangential = std::max(mDamageNormal, mDamageTangential);
        mDamageMoment = std::max(mDamageNormal, mDamageTangential);

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage_parallel_bond::CalculateViscoDamping(double LocalRelVel[3],
                                         double ViscoDampingLocalContactForce[3],
                                         double indentation,
                                         double equiv_visco_damp_coeff_normal,
                                         double equiv_visco_damp_coeff_tangential,
                                         bool& sliding,
                                         int failure_id) {

        KRATOS_TRY

        mUnbondedViscoDampingLocalContactForce[0] = 0.0;
        mUnbondedViscoDampingLocalContactForce[1] = 0.0;
        mUnbondedViscoDampingLocalContactForce[2] = 0.0;
        mBondedViscoDampingLocalContactForce[0] = 0.0;
        mBondedViscoDampingLocalContactForce[1] = 0.0;
        mBondedViscoDampingLocalContactForce[2] = 0.0;

        if (indentation > 0) {
            mUnbondedViscoDampingLocalContactForce[0] = -mUnbondedEquivViscoDampCoeffTangential * LocalRelVel[0];
            mUnbondedViscoDampingLocalContactForce[1] = -mUnbondedEquivViscoDampCoeffTangential * LocalRelVel[1];
            mUnbondedViscoDampingLocalContactForce[2] = -mUnbondedEquivViscoDampCoeffNormal * LocalRelVel[2];
        }

        if (!failure_id) { // Adding bonded and unbonded parts
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

        #ifdef KRATOS_DEBUG
            DemDebugFunctions::CheckIfNan(mUnbondedViscoDampingLocalContactForce, "NAN in Viscous Force in CalculateViscoDamping");
            DemDebugFunctions::CheckIfNan(ViscoDampingLocalContactForce, "NAN in Viscous Force in CalculateViscoDamping");
        #endif

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage_parallel_bond::CalculateViscoDampingCoeff(double& equiv_visco_damp_coeff_normal,
                                              double& equiv_visco_damp_coeff_tangential,
                                              SphericContinuumParticle* element1,
                                              SphericContinuumParticle* element2,
                                              const double kn_el,
                                              const double kt_el) {

        KRATOS_TRY

        const double my_mass    = element1->GetMass();
        const double other_mass = element2->GetMass();

        const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);

        const double my_gamma    = element1->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = element2->GetProperties()[DAMPING_GAMMA];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);

        equiv_visco_damp_coeff_normal     = 2.0 * equiv_gamma * sqrt(equiv_mass * kn_el);
        equiv_visco_damp_coeff_tangential = 2.0 * equiv_gamma * sqrt(equiv_mass * kt_el);

        mUnbondedEquivViscoDampCoeffNormal = 2.0 * equiv_gamma * sqrt(equiv_mass * mUnbondedNormalElasticConstant);
        mUnbondedEquivViscoDampCoeffTangential = 2.0 * equiv_gamma * sqrt(equiv_mass * mUnbondedTangentialElasticConstant);

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage_parallel_bond::ComputeNormalUnbondedForce(double indentation) {

        KRATOS_TRY

        mUnbondedLocalElasticContactForce2 = mUnbondedNormalElasticConstant * indentation;

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage_parallel_bond::CalculateNormalForces(double LocalElasticContactForce[3],
            const double kn_el,
            double equiv_young,
            double indentation,
            double calculation_area,
            double& accumulated_damage,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            int time_steps,
            const ProcessInfo& r_process_info) {

        KRATOS_TRY

        const double tension_limit = 0.5 * (GetContactSigmaMax(element1) + GetContactSigmaMax(element2));
        const double fracture_energy = 0.5 * (element1->GetProperties()[FRACTURE_ENERGY] + element2->GetProperties()[FRACTURE_ENERGY]);
        const double initial_limit_force = tension_limit * calculation_area;

        if (tension_limit) {
            mDamageEnergyCoeff = 2.0 * fracture_energy * kn_el / (calculation_area * tension_limit * tension_limit) - 1.0;
        } else {
            mDamageEnergyCoeff = 0.0;
        }

        KRATOS_ERROR_IF(mDamageEnergyCoeff > 30.0) << "Damage energy is too big!" << std::endl;

        if (mDamageEnergyCoeff < 0.0) {
            mDamageEnergyCoeff = 0.0;
        }

        double k_softening = 0.0;
        double limit_force = 0.0;

        if (mDamageEnergyCoeff) {
            k_softening = kn_el / mDamageEnergyCoeff;
        }

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
        double kn_updated = (1.0 - mDamageNormal) * kn_el;
        double current_normal_force_module = fabs(kn_updated * indentation);
        double delta_accumulated = 0.0;
        if (kn_updated) {
            delta_accumulated = current_normal_force_module / kn_updated;
        }

        double returned_by_mapping_force = 0.0;

        double BondedLocalElasticContactForce2 = 0.0;

        if (indentation >= 0.0) { //COMPRESSION
            if (!failure_type) {
                BondedLocalElasticContactForce2 = kn_updated * indentation;
            } else {
                BondedLocalElasticContactForce2 = 0.0;
            }
        } else { //tension

            if (!failure_type) {
                if (mDamageEnergyCoeff) {
                    limit_force = initial_limit_force * (1.0 + k_softening / kn_el) * kn_updated / (kn_updated + k_softening);
                } else {
                    limit_force = initial_limit_force;
                }

                BondedLocalElasticContactForce2 = kn_updated * indentation;

                if (current_normal_force_module > limit_force) {

                    if (mDamageEnergyCoeff) { // the material can sustain further damage, not failure yet

                        const double delta_at_undamaged_peak = initial_limit_force / kn_el;

                        delta_accumulated = current_normal_force_module / kn_updated;

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

        if (indentation > 0.0) {
           ComputeNormalUnbondedForce(indentation);
        }

        LocalElasticContactForce[2] = BondedLocalElasticContactForce2 + mUnbondedLocalElasticContactForce2;

        double LocalElasticContactTension = LocalElasticContactForce[2] / calculation_area;
        double limit_tension = limit_force / calculation_area;
        double returned_by_mapping_tension = returned_by_mapping_force / calculation_area;
        double current_normal_tension_module = current_normal_force_module / calculation_area;
        double BondedLocalElasticContactTension2 = BondedLocalElasticContactForce2 / calculation_area;
        double initial_limit_tension = initial_limit_force / calculation_area;

        if (mDebugPrintingOption) {

            long unsigned int sphere_id = element1->GetProperties().GetValue(DEBUG_PRINTING_ID_1);
            long unsigned int neigh_sphere_id = element2->GetProperties().GetValue(DEBUG_PRINTING_ID_2);

            if ((element1->Id() == sphere_id) && (element2->Id() == neigh_sphere_id)) {
                std::ofstream normal_forces_file("nl.txt", std::ios_base::out | std::ios_base::app);
                normal_forces_file << r_process_info[TIME] << " " << indentation/*2*/ << " " << LocalElasticContactForce[2]/*3*/ << " "
                                   << limit_force/*4*/ << " " << delta_accumulated/*5*/ << " " << returned_by_mapping_force/*6*/ << " "
                                   << kn_updated/*7*/ << " " << mDamageNormal/*8*/ << " " << failure_type/*9*/ << " "
                                   << current_normal_force_module/*10*/ << " " << mDamageTangential/*11*/ << " "
                                   << BondedLocalElasticContactForce2/*12*/ << " " << mUnbondedLocalElasticContactForce2/*13*/ << " "
                                   << kn_el/*14*/ << " " << mDamageEnergyCoeff/*15*/ << " "
                                   << initial_limit_force/*16*/ << " " << mUnbondedNormalElasticConstant/*17*/ << " "
                                   << LocalElasticContactTension/*18*/ << " " << limit_tension/*19*/ << " "
                                   << returned_by_mapping_tension/*20*/ << " " << current_normal_tension_module/*21*/ << " "
                                   << BondedLocalElasticContactTension2/*22*/ << " " << initial_limit_tension/*23*/ << '\n';
                normal_forces_file.flush();
                normal_forces_file.close();
            }
        }

        #ifdef KRATOS_DEBUG
            DemDebugFunctions::CheckIfNan(BondedLocalElasticContactForce2, "NAN in Bonded Normal Force in CalculateViscoDamping");
            DemDebugFunctions::CheckIfNan(mUnbondedLocalElasticContactForce2, "NAN in Unbonded Normal Force in CalculateViscoDamping");
        #endif

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage_parallel_bond::CalculateTangentialForces(double OldLocalElasticContactForce[3],
            double LocalElasticContactForce[3],
            double LocalElasticExtraContactForce[3],
            double ViscoDampingLocalContactForce[3],
            double LocalCoordSystem[3][3],
            double LocalDeltDisp[3],
            double LocalRelVel[3],
            const double kt_el,
            const double equiv_shear,
            double& contact_sigma,
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

        const double tau_zero = 0.5 * (GetTauZero(element1) + GetTauZero(element2));
        const double internal_friction = 0.5 * (GetInternalFricc(element1) + GetInternalFricc(element2));
        double k_softening = 0.0;
        double tau_strength = 0.0;

        double OldBondedLocalElasticContactForce[2] = {0.0};
        OldBondedLocalElasticContactForce[0] = mBondedScalingFactor * OldLocalElasticContactForce[0];
        OldBondedLocalElasticContactForce[1] = mBondedScalingFactor * OldLocalElasticContactForce[1];

        if (mDamageEnergyCoeff) {
            k_softening = kt_el / mDamageEnergyCoeff;
        }

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
        double kt_updated = (1.0 - mDamageTangential) * kt_el;

        double BondedLocalElasticContactForce[2] = {0.0};

        if (!failure_type) {
            BondedLocalElasticContactForce[0] = OldBondedLocalElasticContactForce[0] - kt_updated * LocalDeltDisp[0]; // 0: first tangential
            BondedLocalElasticContactForce[1] = OldBondedLocalElasticContactForce[1] - kt_updated * LocalDeltDisp[1]; // 1: second tangential

        } else {
            BondedLocalElasticContactForce[0] = 0.0; // 0: first tangential
            BondedLocalElasticContactForce[1] = 0.0; // 1: second tangential
        }

        double delta_accumulated = 0.0;
        double current_tangential_force_module = 0.0;
        double returned_by_mapping_force = 0.0;

        if (!failure_type) { // This means it has not broken yet

            current_tangential_force_module = sqrt(BondedLocalElasticContactForce[0] * BondedLocalElasticContactForce[0]
                                                 + BondedLocalElasticContactForce[1] * BondedLocalElasticContactForce[1]);

            delta_accumulated = current_tangential_force_module / kt_updated;

            returned_by_mapping_force = current_tangential_force_module;

            if (r_process_info[SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION]) { //TODO: use this only for intact bonds (not broken))
                AddContributionOfShearStrainParallelToBond(OldBondedLocalElasticContactForce, LocalElasticExtraContactForce, element1->mNeighbourElasticExtraContactForces[i_neighbour_count], LocalCoordSystem, kt_el, calculation_area,  element1, element2);
            }

            contact_sigma = LocalElasticContactForce[2] / calculation_area;
            contact_tau = current_tangential_force_module / calculation_area;

            double updated_max_tau_strength = tau_zero;

            if (contact_sigma >= 0) {
                updated_max_tau_strength += internal_friction * contact_sigma;
            }
            tau_strength = updated_max_tau_strength * (1.0 + k_softening / kt_el) * kt_updated / (kt_updated + k_softening);

            if (contact_tau > tau_strength) { // damage

                if (mDamageEnergyCoeff) { // the material can sustain further damage, not failure yet

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

        // HERE IT STARTS THE UNBONDED PART
        double OldUnbondedLocalElasticContactForce[2] = {0.0};
        double UnbondedLocalElasticContactForce[2] = {0.0};
        double ActualTotalShearForce = 0.0;
        double max_admissible_shear_force = 0.0;
        double fraction = 0.0;

        if (indentation <= 0.0) {
            UnbondedLocalElasticContactForce[0] = 0.0;
            UnbondedLocalElasticContactForce[1] = 0.0;
        } else {
            OldUnbondedLocalElasticContactForce[0] = mUnbondedScalingFactor * OldLocalElasticContactForce[0];
            OldUnbondedLocalElasticContactForce[1] = mUnbondedScalingFactor * OldLocalElasticContactForce[1];

            UnbondedLocalElasticContactForce[0] = OldUnbondedLocalElasticContactForce[0] - mUnbondedTangentialElasticConstant * LocalDeltDisp[0];
            UnbondedLocalElasticContactForce[1] = OldUnbondedLocalElasticContactForce[1] - mUnbondedTangentialElasticConstant * LocalDeltDisp[1];


            const double my_tg_of_static_friction_angle        = element1->GetTgOfStaticFrictionAngle();
            const double neighbour_tg_of_static_friction_angle = element2->GetTgOfStaticFrictionAngle();
            const double equiv_tg_of_static_fri_ang            = 0.5 * (my_tg_of_static_friction_angle + neighbour_tg_of_static_friction_angle);

            const double my_tg_of_dynamic_friction_angle        = element1->GetTgOfDynamicFrictionAngle();
            const double neighbour_tg_of_dynamic_friction_angle = element2->GetTgOfDynamicFrictionAngle();
            const double equiv_tg_of_dynamic_fri_ang            = 0.5 * (my_tg_of_dynamic_friction_angle + neighbour_tg_of_dynamic_friction_angle);
            const double my_friction_decay_coefficient          = element1->GetFrictionDecayCoefficient();
            const double neighbour_friction_decay_coefficient   = element2->GetProperties()[FRICTION_DECAY];
            double equiv_friction_decay_coefficient             = 0.5 * (my_friction_decay_coefficient + neighbour_friction_decay_coefficient);

            if(equiv_tg_of_static_fri_ang < 0.0 || equiv_tg_of_dynamic_fri_ang < 0.0) {
                KRATOS_ERROR << "The averaged friction is negative for one contact of element with Id: "<< element1->Id()<<std::endl;
            }

            const double ShearRelVel = sqrt(LocalRelVel[0] * LocalRelVel[0] + LocalRelVel[1] * LocalRelVel[1]);
            double equiv_friction = equiv_tg_of_dynamic_fri_ang + (equiv_tg_of_static_fri_ang - equiv_tg_of_dynamic_fri_ang) * exp(-equiv_friction_decay_coefficient * ShearRelVel);

            double max_admissible_shear_force = (mUnbondedLocalElasticContactForce2 + mUnbondedViscoDampingLocalContactForce[2]) * equiv_friction;

            if (equiv_tg_of_static_fri_ang < 0.0 || equiv_tg_of_dynamic_fri_ang < 0.0) {
                KRATOS_ERROR << "The averaged friction is negative for one contact of element with Id: "<< element1->Id()<<std::endl;
            }

            const double tangential_contact_force_0 = UnbondedLocalElasticContactForce[0] + mUnbondedViscoDampingLocalContactForce[0];
            const double tangential_contact_force_1 = UnbondedLocalElasticContactForce[1] + mUnbondedViscoDampingLocalContactForce[1];

            ActualTotalShearForce = sqrt(tangential_contact_force_0 * tangential_contact_force_0 + tangential_contact_force_1 * tangential_contact_force_1);

            if (ActualTotalShearForce > max_admissible_shear_force) {
                const double ActualElasticShearForce = sqrt(UnbondedLocalElasticContactForce[0] * UnbondedLocalElasticContactForce[0] + UnbondedLocalElasticContactForce[1] * UnbondedLocalElasticContactForce[1]);

                const double dot_product = UnbondedLocalElasticContactForce[0] * mUnbondedViscoDampingLocalContactForce[0] + UnbondedLocalElasticContactForce[1] * mUnbondedViscoDampingLocalContactForce[1];
                const double ViscoDampingLocalContactForceModule = sqrt(mUnbondedViscoDampingLocalContactForce[0] * mUnbondedViscoDampingLocalContactForce[0] +
                                                                        mUnbondedViscoDampingLocalContactForce[1] * mUnbondedViscoDampingLocalContactForce[1]);

                if (dot_product >= 0.0) {
                    if (ActualElasticShearForce > max_admissible_shear_force) {
                        fraction = max_admissible_shear_force / ActualElasticShearForce;
                        UnbondedLocalElasticContactForce[0]      = UnbondedLocalElasticContactForce[0] * fraction;
                        UnbondedLocalElasticContactForce[1]      = UnbondedLocalElasticContactForce[1] * fraction;
                        mUnbondedViscoDampingLocalContactForce[0] = 0.0;
                        mUnbondedViscoDampingLocalContactForce[1] = 0.0;
                    }
                    else {
                        const double ActualViscousShearForce = max_admissible_shear_force - ActualElasticShearForce;
                        fraction = ActualViscousShearForce / ViscoDampingLocalContactForceModule;
                        mUnbondedViscoDampingLocalContactForce[0] *= fraction;
                        mUnbondedViscoDampingLocalContactForce[1] *= fraction;
                    }
                }
                else {
                    if (ViscoDampingLocalContactForceModule >= ActualElasticShearForce) {
                        fraction = (max_admissible_shear_force + ActualElasticShearForce) / ViscoDampingLocalContactForceModule;
                        mUnbondedViscoDampingLocalContactForce[0] *= fraction;
                        mUnbondedViscoDampingLocalContactForce[1] *= fraction;
                    }
                    else {
                        fraction = max_admissible_shear_force / ActualElasticShearForce;
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

        double local_elastic_force_modulus = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] +
                                                  LocalElasticContactForce[1] * LocalElasticContactForce[1]);

        if (local_elastic_force_modulus) {
            mBondedScalingFactor = (BondedLocalElasticContactForce[0] * LocalElasticContactForce[0] +
                                    BondedLocalElasticContactForce[1] * LocalElasticContactForce[1]) / (local_elastic_force_modulus * local_elastic_force_modulus);
            mUnbondedScalingFactor = (UnbondedLocalElasticContactForce[0] * LocalElasticContactForce[0] +
                                      UnbondedLocalElasticContactForce[1] * LocalElasticContactForce[1]) / (local_elastic_force_modulus * local_elastic_force_modulus);
        } else {
            mBondedScalingFactor = mUnbondedScalingFactor = 0.0;
        }
        double returned_by_mapping_tension = returned_by_mapping_force / calculation_area;

        if (mDebugPrintingOption) {
            double local_elastic_force_modulus_bonded_only = sqrt(BondedLocalElasticContactForce[0] * BondedLocalElasticContactForce[0] +
                                                              BondedLocalElasticContactForce[1] * BondedLocalElasticContactForce[1]);
            double local_elastic_force_modulus_unbonded_only = sqrt(UnbondedLocalElasticContactForce[0] * UnbondedLocalElasticContactForce[0] +
                                                              UnbondedLocalElasticContactForce[1] * UnbondedLocalElasticContactForce[1]);
            long unsigned int sphere_id = element1->GetProperties().GetValue(DEBUG_PRINTING_ID_1);
            long unsigned int neigh_sphere_id = element2->GetProperties().GetValue(DEBUG_PRINTING_ID_2);
            double quotient = local_elastic_force_modulus / calculation_area;
            double quotient_bonded_only = local_elastic_force_modulus_bonded_only / calculation_area;
            double quotient_unbonded_only = local_elastic_force_modulus_unbonded_only / calculation_area;

            if ((element1->Id() == sphere_id) && (element2->Id() == neigh_sphere_id)) {
                std::ofstream tangential_forces_file("tg.txt", std::ios_base::out | std::ios_base::app);
                tangential_forces_file << r_process_info[TIME] << " " << int(failure_type)/*2*/ << " " << tau_strength/*3*/ << " "
                                    << kt_updated/*4*/ << " " << int(sliding)/*5*/ << " "
                                    << contact_sigma/*6*/ << " " << mDamageNormal/*7*/ << " " << contact_tau/*8*/ << " "
                                    << max_admissible_shear_force/*9*/ << " " << mDamageTangential/*10*/ << " " << LocalElasticContactForce[0]/*11*/ << " "
                                    << OldBondedLocalElasticContactForce[0]/*12*/ << " " << mUnbondedLocalElasticContactForce2/*13*/ << " "
                                    << BondedLocalElasticContactForce[0]/*14*/ << " " << mBondedScalingFactor/*15*/ << " "
                                    << UnbondedLocalElasticContactForce[0]/*16*/ << " " << LocalDeltDisp[0]/*17*/ << " " << mUnbondedScalingFactor/*18*/ << " "
                                    << OldLocalElasticContactForce[0]/*19*/ << " " << LocalDeltDisp[1]/*20*/ << " " << LocalElasticContactForce[1]/*21*/ << " "
                                    << BondedLocalElasticContactForce[1]/*22*/ << " " << UnbondedLocalElasticContactForce[1]/*23*/ << " "
                                    << OldBondedLocalElasticContactForce[1]/*24*/ << " " << OldUnbondedLocalElasticContactForce[1]/*25*/ << " "
                                    << delta_accumulated/*26*/ << " " << current_tangential_force_module/*27*/ << " "
                                    << returned_by_mapping_force/*28*/ << " " << quotient/*29*/ << " "
                                    << quotient_bonded_only/*30*/ << " " << quotient_unbonded_only/*31*/ << " "
                                    << returned_by_mapping_tension/*32*/ << " " << tau_strength * calculation_area /*33*/ << " "
                                    << mDamageEnergyCoeff/*34*/ << " " << LocalElasticContactForce[2]/*35*/ << " "
                                    << indentation/*36*/<< '\n';
                tangential_forces_file.flush();
                tangential_forces_file.close();
            }
        }

        KRATOS_CATCH("")
    }

    double DEM_KDEM_with_damage_parallel_bond::LocalMaxSearchDistance(const int i,
                                            SphericContinuumParticle* element1,
                                            SphericContinuumParticle* element2) {

        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();
        double tension_limit;

        // calculation of equivalent Young modulus
        double myYoung = element1->GetProperties()[BONDED_MATERIAL_YOUNG_MODULUS];
        double other_young = element2->GetProperties()[BONDED_MATERIAL_YOUNG_MODULUS];
        double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);

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

        if (&element1_props == &element2_props) {
            tension_limit = GetContactSigmaMax(element1);
        } else {
            tension_limit = 0.5 * (GetContactSigmaMax(element1) + GetContactSigmaMax(element2));
        }

        const double Ntstr_el = tension_limit * calculation_area;
        double u1 = Ntstr_el / kn_el;
        if (u1 > 2.0 * radius_sum) {
            u1 = 2.0 * radius_sum;
        } // avoid error in special cases with too high tensile
        return u1;
    }

} // namespace Kratos
