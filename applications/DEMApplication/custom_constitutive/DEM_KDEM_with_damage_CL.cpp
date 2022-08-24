// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_with_damage_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_with_damage::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_with_damage(*this));
        return p_clone;
    }

    void DEM_KDEM_with_damage::Check(Properties::Pointer pProp) const {

        if (!pProp->Has(SHEAR_ENERGY_COEF)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable SHEAR_ENERGY_COEF should be present in the properties when using DEM_KDEM_with_damage. A default value of 0.0 was assigned."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(SHEAR_ENERGY_COEF) = 0.0;
        }
    }

    void DEM_KDEM_with_damage::CalculateForces(const ProcessInfo& r_process_info,
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
        
        CalculateNormalAndTangentialDamageComponents();

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage::CalculateNormalForces(double LocalElasticContactForce[3],
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

        const double tension_limit = GetContactSigmaMax();
        const double& damage_energy_coeff = (*mpProperties)[SHEAR_ENERGY_COEF];
        double k_unload = 0.0;
        double limit_force = 0.0;

        if (damage_energy_coeff) {
            k_unload = kn_el / damage_energy_coeff;
        }

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
        double kn_updated = (1.0 - mDamageNormal) * kn_el;
        double current_normal_force_module = fabs(kn_updated * indentation);
        double delta_acummulated = current_normal_force_module / kn_updated;
        double returned_by_mapping_force = current_normal_force_module;

        if (indentation >= 0.0) { //COMPRESSION
            LocalElasticContactForce[2] = kn_el * indentation;
        } else { //tension
            if (!failure_type) {

                const double initial_limit_force = tension_limit * calculation_area;
                limit_force = (1.0 - mDamageNormal) * initial_limit_force;
                LocalElasticContactForce[2] = kn_updated * indentation;

                if ((current_normal_force_module > limit_force) && !(*mpProperties)[IS_UNBREAKABLE]) {

                    if (!damage_energy_coeff) { // there is no damage energy left
                        failure_type = 4; // failure by traction
                    } else { // the material can sustain further damage, not failure yet
                        const double delta_at_undamaged_peak = initial_limit_force / kn_el;

                        if (kn_updated) {
                            delta_acummulated = current_normal_force_module / kn_updated;
                        } else {
                            delta_acummulated = delta_at_undamaged_peak + initial_limit_force / k_unload;
                        }

                        returned_by_mapping_force = initial_limit_force - k_unload * (delta_acummulated - delta_at_undamaged_peak);

                        if (returned_by_mapping_force < 0.0) {
                            returned_by_mapping_force = 0.0;
                        }

                        LocalElasticContactForce[2] = -returned_by_mapping_force;

                        mDamageNormal = 1.0 - (returned_by_mapping_force / delta_acummulated) / kn_el;
                        if (mDamageNormal > mDamageThresholdTolerance) {
                            failure_type = 4; // failure by traction
                        }
                    }
                }
            }
            else {
                LocalElasticContactForce[2] = 0.0;
            }
        }

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage::CalculateTangentialForces(double OldLocalElasticContactForce[3],
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

        const double& tau_zero = (*mpProperties)[CONTACT_TAU_ZERO];
        const double& internal_friction = (*mpProperties)[CONTACT_INTERNAL_FRICC];
        const double& damage_energy_coeff = (*mpProperties)[SHEAR_ENERGY_COEF];
        double k_unload = 0.0;
        double tau_strength = 0.0;

        if (damage_energy_coeff) {
            k_unload = kt_el / damage_energy_coeff;
        }

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
        double kt_updated = (1.0 - mDamageTangential) * kt_el;

        if (!failure_type) {
            LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - kt_updated * LocalDeltDisp[0];
            LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - kt_updated * LocalDeltDisp[1];
        } else {
            LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - kt_el * LocalDeltDisp[0];
            LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - kt_el * LocalDeltDisp[1];
        }

        double current_tangential_force_module = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                                    + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

        double delta_acummulated = current_tangential_force_module / kt_updated;

        double returned_by_mapping_force = current_tangential_force_module;

        if (!failure_type) { // This means it has not broken yet

            if (r_process_info[SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION]) { //TODO: use this only for intact bonds (not broken))
                AddContributionOfShearStrainParallelToBond(OldLocalElasticContactForce, LocalElasticExtraContactForce, element1->mNeighbourElasticExtraContactForces[i_neighbour_count], LocalCoordSystem, kt_el, calculation_area,  element1, element2);
            }

            contact_sigma = LocalElasticContactForce[2] / calculation_area;
            contact_tau = current_tangential_force_module / calculation_area;

            double updated_max_tau_strength = tau_zero;
            tau_strength = (1.0 - mDamageTangential) * tau_zero;

            if (contact_sigma >= 0) {
                tau_strength += (1.0 - mDamageTangential) * internal_friction * contact_sigma;
                updated_max_tau_strength += internal_friction * contact_sigma;
            }

            if ((contact_tau > tau_strength) && !(*mpProperties)[IS_UNBREAKABLE]) { // damage

                if (!damage_energy_coeff) { // there is no damage energy left
                    failure_type = 2; // failure by shear
                } else { // the material can sustain further damage, not failure yet
                    const double delta_at_undamaged_peak = updated_max_tau_strength * calculation_area / kt_el;

                    if (kt_updated) {
                        delta_acummulated = current_tangential_force_module / kt_updated;
                    } else {
                        delta_acummulated = delta_at_undamaged_peak + updated_max_tau_strength * calculation_area / k_unload;
                    }

                    returned_by_mapping_force = updated_max_tau_strength * calculation_area - k_unload * (delta_acummulated - delta_at_undamaged_peak);

                    if (returned_by_mapping_force < 0.0) {
                        returned_by_mapping_force = 0.0;
                    }

                    if (current_tangential_force_module) {
                        LocalElasticContactForce[0] = (returned_by_mapping_force / current_tangential_force_module) * LocalElasticContactForce[0];
                        LocalElasticContactForce[1] = (returned_by_mapping_force / current_tangential_force_module) * LocalElasticContactForce[1];
                    }

                    current_tangential_force_module = returned_by_mapping_force; // computed only for printing purposes

                    mDamageTangential = 1.0 - (returned_by_mapping_force / delta_acummulated) / kt_el;

                    if (mDamageTangential > mDamageThresholdTolerance) {
                        failure_type = 2; // failure by shear
                    }
                }
            }
        } else {
            const double& equiv_tg_of_static_fri_ang = (*mpProperties)[STATIC_FRICTION];
            const double& equiv_tg_of_dynamic_fri_ang = (*mpProperties)[DYNAMIC_FRICTION];
            const double& equiv_friction_decay_coefficient = (*mpProperties)[FRICTION_DECAY];

            const double ShearRelVel = sqrt(LocalRelVel[0] * LocalRelVel[0] + LocalRelVel[1] * LocalRelVel[1]);
            double equiv_friction = equiv_tg_of_dynamic_fri_ang + (equiv_tg_of_static_fri_ang - equiv_tg_of_dynamic_fri_ang) * exp(-equiv_friction_decay_coefficient * ShearRelVel);

            //
            double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];

            if (normal_contact_force < 0.0) {
                normal_contact_force = 0.0;
                ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
            }

            double maximum_admissible_shear_force = normal_contact_force * equiv_friction;

            const double tangential_contact_force_0 = LocalElasticContactForce[0] + ViscoDampingLocalContactForce[0];
            const double tangential_contact_force_1 = LocalElasticContactForce[1] + ViscoDampingLocalContactForce[1];

            const double ActualTotalShearForce = sqrt(tangential_contact_force_0 * tangential_contact_force_0 + tangential_contact_force_1 * tangential_contact_force_1);

            if (ActualTotalShearForce > maximum_admissible_shear_force) {

                const double ActualElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

                const double dot_product = LocalElasticContactForce[0] * ViscoDampingLocalContactForce[0] + LocalElasticContactForce[1] * ViscoDampingLocalContactForce[1];
                const double ViscoDampingLocalContactForceModule = sqrt(ViscoDampingLocalContactForce[0] * ViscoDampingLocalContactForce[0] +\
                                                                        ViscoDampingLocalContactForce[1] * ViscoDampingLocalContactForce[1]);

                if (dot_product >= 0.0) {

                    if (ActualElasticShearForce > maximum_admissible_shear_force) {
                        const double fraction = maximum_admissible_shear_force / ActualElasticShearForce;
                        LocalElasticContactForce[0]      = LocalElasticContactForce[0] * fraction;
                        LocalElasticContactForce[1]      = LocalElasticContactForce[1] * fraction;
                        ViscoDampingLocalContactForce[0] = 0.0;
                        ViscoDampingLocalContactForce[1] = 0.0;
                    }
                    else {
                        const double ActualViscousShearForce = maximum_admissible_shear_force - ActualElasticShearForce;
                        const double fraction = ActualViscousShearForce / ViscoDampingLocalContactForceModule;
                        ViscoDampingLocalContactForce[0] *= fraction;
                        ViscoDampingLocalContactForce[1] *= fraction;
                    }
                }
                else {
                    if (ViscoDampingLocalContactForceModule >= ActualElasticShearForce) {
                        const double fraction = (maximum_admissible_shear_force + ActualElasticShearForce) / ViscoDampingLocalContactForceModule;
                        ViscoDampingLocalContactForce[0] *= fraction;
                        ViscoDampingLocalContactForce[1] *= fraction;
                    }
                    else {
                        const double fraction = maximum_admissible_shear_force / ActualElasticShearForce;
                        LocalElasticContactForce[0]      = LocalElasticContactForce[0] * fraction;
                        LocalElasticContactForce[1]      = LocalElasticContactForce[1] * fraction;
                        ViscoDampingLocalContactForce[0] = 0.0;
                        ViscoDampingLocalContactForce[1] = 0.0;
                    }
                }
                sliding = true;
            }
        }

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage::CalculateNormalAndTangentialDamageComponents() {

        KRATOS_TRY

        mDamageReal += std::sqrt((mDamageNormal - mDamageReal) * (mDamageNormal - mDamageReal) + (mDamageTangential - mDamageReal) * (mDamageTangential - mDamageReal));
        mDamageNormal = mDamageReal;
        mDamageTangential = mDamageReal;
        mDamageMoment = mDamageReal;

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
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

        DEM_KDEM_soft_torque::ComputeParticleRotationalMoments(element,
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

        ElasticLocalRotationalMoment[0] *= (1.0 - mDamageMoment);
        ElasticLocalRotationalMoment[1] *= (1.0 - mDamageMoment);
        ElasticLocalRotationalMoment[2] *= (1.0 - mDamageMoment);

        KRATOS_CATCH("")
    }//ComputeParticleRotationalMoments

} // namespace Kratos
