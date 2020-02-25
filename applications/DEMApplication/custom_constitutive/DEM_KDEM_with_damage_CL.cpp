// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_with_damage_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    void DEM_KDEM_with_damage::Initialize(SphericContinuumParticle* element) {}

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_with_damage::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_with_damage(*this));
        return p_clone;
    }

    void DEM_KDEM_with_damage::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_with_damage to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
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
                                int search_control,
                                DenseVector<int>& search_control_vector,
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

        CalculateTangentialForces(OldLocalElasticContactForce,
                LocalElasticContactForce,
                LocalElasticExtraContactForce,
                LocalCoordSystem,
                LocalDeltDisp,
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
                search_control,
                search_control_vector,
                r_process_info);

        FindMaximumValueOfNormalAndTangentialDamageComponents();

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

        const double tension_limit = 0.5 * (GetContactSigmaMax(element1) + GetContactSigmaMax(element2));
        const double damage_energy_coeff = 0.5 * (element1->GetProperties()[SHEAR_ENERGY_COEF] + element2->GetProperties()[SHEAR_ENERGY_COEF]);
        double k_unload = 0.0;
        double limit_force = 0.0;
        static bool first_time_entered = true;
        const unsigned int sphere_id = 22222222;
        const std::string filename = "normal_forces_damage.txt";

        if (element1->Id() == sphere_id) {
            static std::ifstream ifile(filename.c_str());
            if ((bool) ifile && first_time_entered) {
                std::remove("normal_forces_damage.txt");
                first_time_entered = false;
            }
        }

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

                if (current_normal_force_module > limit_force) {

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

        if (element1->Id() == sphere_id) {
            static std::ofstream normal_forces_file("normal_forces_damage.txt", std::ios_base::out | std::ios_base::app);
            normal_forces_file << r_process_info[TIME] << " " << indentation << " " << LocalElasticContactForce[2] << " " << limit_force << " "
                               << delta_acummulated << " " << returned_by_mapping_force << " " << kn_updated << " " << mDamageNormal << " "
                               << failure_type << " " << current_normal_force_module << " " << mDamageTangential <<'\n';
            normal_forces_file.flush();
        }

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage::CalculateTangentialForces(double OldLocalElasticContactForce[3],
            double LocalElasticContactForce[3],
            double LocalElasticExtraContactForce[3],
            double LocalCoordSystem[3][3],
            double LocalDeltDisp[3],
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
            int search_control,
            DenseVector<int>& search_control_vector,
            const ProcessInfo& r_process_info) {

        KRATOS_TRY

        const double tau_zero = 0.5 * (GetTauZero(element1) + GetTauZero(element2));
        const double internal_friction = 0.5 * (GetInternalFricc(element1) + GetInternalFricc(element2));
        const double damage_energy_coeff = 0.5 * (element1->GetProperties()[SHEAR_ENERGY_COEF] + element2->GetProperties()[SHEAR_ENERGY_COEF]);
        double k_unload = 0.0;
        double tau_strength = 0.0;
        static bool first_time_entered = true;
        int damage_process = 0;
        double equiv_tg_of_fri_ang;
        double maximum_frictional_shear_force = 0.0;
        const unsigned int sphere_id = 22222222;
        const std::string filename = "tangential_forces_damage.txt";

        if (element1->Id() == sphere_id) {
            static std::ifstream ifile(filename.c_str());
            if ((bool) ifile && first_time_entered) {
                std::remove("tangential_forces_damage.txt");
                first_time_entered = false;
            }
        }

        if (damage_energy_coeff) {
            k_unload = kt_el / damage_energy_coeff;
        }

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
        double kt_updated = (1.0 - mDamageTangential) * kt_el;

        if (!failure_type) {
            LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - kt_updated * LocalDeltDisp[0]; // 0: first tangential
            LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - kt_updated * LocalDeltDisp[1]; // 1: second tangential
        } else {
            LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - kt_el * LocalDeltDisp[0]; // 0: first tangential
            LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - kt_el * LocalDeltDisp[1]; // 1: second tangential
        }
        double total_delta_displ_module = sqrt(LocalDeltDisp[0] * LocalDeltDisp[0] + LocalDeltDisp[1] * LocalDeltDisp[1]);

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

            if (contact_tau > tau_strength) { // damage

                damage_process = 1;

                if (!damage_energy_coeff) { // there is no damage energy left
                    failure_type = 2; // failure by shear
                } else { // the material can sustain further damage, not failure yet
                    const double delta_at_undamaged_peak = updated_max_tau_strength * calculation_area / kt_el;

                    if (kt_updated) {
                        delta_acummulated = current_tangential_force_module / kt_updated;
                    } else {
                        delta_acummulated = delta_at_undamaged_peak + updated_max_tau_strength * calculation_area / k_unload;
                    }
                    if (element1->Id() == sphere_id) {
                        KRATOS_WATCH(LocalDeltDisp[0])
                        KRATOS_WATCH(OldLocalElasticContactForce[0]/kt_updated)
                        KRATOS_WATCH((LocalElasticContactForce[0]-OldLocalElasticContactForce[0])/kt_updated)
                        KRATOS_WATCH(LocalElasticContactForce[0])
                        KRATOS_WATCH(kt_updated)
                        KRATOS_WATCH(delta_at_undamaged_peak)
                        KRATOS_WATCH(delta_acummulated)
                    }

                    returned_by_mapping_force = updated_max_tau_strength * calculation_area - k_unload * (delta_acummulated - delta_at_undamaged_peak);

                    if (element1->Id() == sphere_id) {
                        KRATOS_WATCH(delta_acummulated - delta_at_undamaged_peak)
                        KRATOS_WATCH(returned_by_mapping_force)
                    }

                    if (returned_by_mapping_force < 0.0) {
                        returned_by_mapping_force = 0.0;
                    }

                    if (current_tangential_force_module) {
                        LocalElasticContactForce[0] = (returned_by_mapping_force / current_tangential_force_module) * LocalElasticContactForce[0];
                        LocalElasticContactForce[1] = (returned_by_mapping_force / current_tangential_force_module) * LocalElasticContactForce[1];
                    }

                    current_tangential_force_module = returned_by_mapping_force; // computed only for printing purposes

                    mDamageTangential = 1.0 - (returned_by_mapping_force / delta_acummulated) / kt_el;
                    if (element1->Id() == sphere_id) {
                        KRATOS_WATCH(mDamageTangential)
                    }

                    if (mDamageTangential > mDamageThresholdTolerance) {
                        failure_type = 2; // failure by shear
                    }
                }
            }
        } else {
            equiv_tg_of_fri_ang = 0.5 * (element1->GetTgOfFrictionAngle() + element2->GetTgOfFrictionAngle());

            if(equiv_tg_of_fri_ang < 0.0) {
                KRATOS_ERROR << "The averaged friction is negative for one contact of element with Id: "<< element1->Id()<<std::endl;
            }

            maximum_frictional_shear_force = equiv_tg_of_fri_ang * LocalElasticContactForce[2];

            if (maximum_frictional_shear_force < 0.0) maximum_frictional_shear_force = 0.0;

            if (current_tangential_force_module > maximum_frictional_shear_force) {
                LocalElasticContactForce[0] = (maximum_frictional_shear_force / current_tangential_force_module) * LocalElasticContactForce[0];
                LocalElasticContactForce[1] = (maximum_frictional_shear_force / current_tangential_force_module) * LocalElasticContactForce[1];
                sliding = true;
            }
        }

        if (element1->Id() == sphere_id) {
            static std::ofstream tangential_forces_file("tangential_forces_damage.txt", std::ios_base::out | std::ios_base::app);
            tangential_forces_file << r_process_info[TIME] << " " << int(failure_type) << " " << LocalElasticContactForce[0] << " "
                                   << tau_strength << " " << delta_acummulated << " " << returned_by_mapping_force << " "
                                   << kt_updated << " " << damage_process << " " << int(sliding) << " " << contact_sigma << " " << mDamageNormal << " "
                                   << contact_tau << " " << current_tangential_force_module << " " << LocalElasticContactForce[2] << " "
                                   << maximum_frictional_shear_force << " " << mDamageTangential << " " << LocalDeltDisp[0] << " "
                                   << total_delta_displ_module << '\n';
            tangential_forces_file.flush();
        }

        KRATOS_CATCH("")
    }

    void DEM_KDEM_with_damage::FindMaximumValueOfNormalAndTangentialDamageComponents() {

        KRATOS_TRY

        mDamageNormal = std::max(mDamageNormal, mDamageTangential);
        mDamageTangential = std::max(mDamageNormal, mDamageTangential);

        KRATOS_CATCH("")
    }

} // namespace Kratos
