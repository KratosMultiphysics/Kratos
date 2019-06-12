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

    void DEM_KDEM_with_damage::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
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

        const double tension_limit = 0.5 * 1e6 * (element1->GetFastProperties()->GetContactSigmaMin() + element2->GetFastProperties()->GetContactSigmaMin());
        const double damage_energy = 0.5 * (element1->GetProperties()[SHEAR_ENERGY_COEF] + element2->GetProperties()[SHEAR_ENERGY_COEF]);
        const double damage_threshold_tolerance = 0.9;
        double k_unload = 0.0;
        double limit_force = 0.0;
        double delta_acummulated = 0.0;
        double returned_by_mapping_force = 0.0;
        static bool first_time_entered = true;

        const std::string filename = "normal_forces_damage.txt";
        std::ifstream ifile(filename.c_str());

        if (element1->Id() == 222222222222) {

            if ((bool) ifile && first_time_entered) {
                //std::remove("normal_forces_damage.txt");
                first_time_entered = false;
            }
        }

        if (damage_energy) {
            k_unload = kn_el / damage_energy;
        }

        if (mKnUpdated > kn_el) { // This only happens the first time we enter the function
            mKnUpdated = kn_el;
        }

        if (indentation >= 0.0) { //COMPRESSION
            LocalElasticContactForce[2] = kn_el * indentation;
        } else { //tension

            delta_acummulated = -indentation;

            int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

            if (!failure_type) {

                const double initial_limit_force = tension_limit * calculation_area;
                limit_force = (1.0 - mDamage) * tension_limit * calculation_area;
                LocalElasticContactForce[2] = mKnUpdated * indentation;

                if (fabs(LocalElasticContactForce[2]) > limit_force) {

                    if (!damage_energy) {
                        failure_type = 4; // failure by traction
                        LocalElasticContactForce[2] = 0.0;
                    } else {
                        const double delta_at_undamaged_peak = initial_limit_force / kn_el;

                        if (mKnUpdated) {
                            delta_acummulated = fabs(LocalElasticContactForce[2]) / mKnUpdated;
                        } else {
                            delta_acummulated = delta_at_undamaged_peak + initial_limit_force / k_unload;
                        }

                        returned_by_mapping_force = initial_limit_force - k_unload * (delta_acummulated - delta_at_undamaged_peak);

                        if (returned_by_mapping_force < 0.0) {
                            returned_by_mapping_force = 0.0;
                        }

                        LocalElasticContactForce[2] = -returned_by_mapping_force;

                        if (delta_acummulated) {
                            mKnUpdated = returned_by_mapping_force / delta_acummulated;
                        }

                        mDamage = 1.0 - returned_by_mapping_force/initial_limit_force;
                        if (mDamage > damage_threshold_tolerance) {
                            failure_type = 4; // failure by traction
                        }
                    }
                }
            }
            else {
                LocalElasticContactForce[2] = 0.0;
            }
        }

        if (element1->Id() == 222222222222) {
            static std::ofstream normal_forces_file("normal_forces_damage.txt", std::ios_base::out | std::ios_base::app);
            normal_forces_file << r_process_info[TIME] << " " << indentation << " " << LocalElasticContactForce[2] << " " << limit_force << " "
                               << delta_acummulated << " " << returned_by_mapping_force << " " << mKnUpdated << '\n';
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

        const double mTauZero = 0.5 * 1e6 * (GetTauZero(element1) + GetTauZero(element2));
        const double mInternalFriction = 0.5 * (GetInternalFricc(element1) + GetInternalFricc(element2));
        const double damage_energy = 0.5 * (element1->GetProperties()[SHEAR_ENERGY_COEF] + element2->GetProperties()[SHEAR_ENERGY_COEF]);
        const double damage_threshold_tolerance = 0.9;
        double k_unload = 0.0;
        double tau_strength = 0.0;
        static bool first_time_entered = true;
        int failure = 0;
        int damage = 0;
        int slide = 0;
        double equiv_tg_of_fri_ang;
        double maximum_frictional_shear_force;

        const std::string filename = "tangential_forces_damage.txt";
        std::ifstream ifile(filename.c_str());

        if (element1->Id() == 222222222222) {

            if ((bool) ifile && first_time_entered) {
                //std::remove("tangential_forces_damage.txt");
                first_time_entered = false;
            }
        }

        if (damage_energy) {
            k_unload = kt_el / damage_energy;
        }

        if (mKtUpdated > kt_el) { // This only happens the first time we enter the function
            mKtUpdated = kt_el;
        }

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - mKtUpdated * LocalDeltDisp[0]; // 0: first tangential
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - mKtUpdated * LocalDeltDisp[1]; // 1: second tangential

        const double current_tangential_force_module = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                                          + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

        double returned_by_mapping_force = current_tangential_force_module;
        double delta_acummulated = returned_by_mapping_force / kt_el;

        if (element1->Id() == 222222222222) {
            KRATOS_WATCH(r_process_info[TIME])
            KRATOS_WATCH(kt_el)
            KRATOS_WATCH(returned_by_mapping_force)
            KRATOS_WATCH(delta_acummulated)
        }

        if (!failure_type) { // This means it has not broken

            if (r_process_info[SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION]) { //TODO: use this only for intact bonds (not broken))
                AddContributionOfShearStrainParallelToBond(OldLocalElasticContactForce, LocalElasticExtraContactForce, element1->mNeighbourElasticExtraContactForces[i_neighbour_count], LocalCoordSystem, kt_el, calculation_area,  element1, element2);
            }

            contact_sigma = LocalElasticContactForce[2] / calculation_area;
            contact_tau = current_tangential_force_module / calculation_area;

            const double initial_tau_strength = mTauZero;
            tau_strength = (1.0 - mDamage) * mTauZero;

            if (contact_sigma >= 0) {
                tau_strength += (1.0 - mDamage) * mInternalFriction * contact_sigma;
            }

            if (contact_tau > tau_strength) { // damage

                damage = 1;

                if (!damage_energy) {
                    failure_type = 2; // failure by shear
                    //
                    equiv_tg_of_fri_ang = 0.5 * (element1->GetTgOfFrictionAngle() + element2->GetTgOfFrictionAngle());
                    maximum_frictional_shear_force = equiv_tg_of_fri_ang * LocalElasticContactForce[2];

                    if (maximum_frictional_shear_force < 0.0) maximum_frictional_shear_force = 0.0;

                    if (current_tangential_force_module > maximum_frictional_shear_force) {
                        LocalElasticContactForce[0] = (maximum_frictional_shear_force / current_tangential_force_module) * LocalElasticContactForce[0];
                        LocalElasticContactForce[1] = (maximum_frictional_shear_force / current_tangential_force_module) * LocalElasticContactForce[1];
                        sliding = true;
                        slide = 1;
                    }
                    //
                } else {
                    const double delta_at_undamaged_peak = initial_tau_strength * calculation_area / kt_el;

                    if (mKtUpdated) {
                        delta_acummulated = current_tangential_force_module / mKtUpdated;
                    } else {
                        delta_acummulated = delta_at_undamaged_peak + initial_tau_strength * calculation_area / k_unload;
                    }
                    if (element1->Id() == 222222222222) {
                        KRATOS_WATCH(r_process_info[TIME])
                        KRATOS_WATCH(current_tangential_force_module)
                        KRATOS_WATCH(mKtUpdated)
                        KRATOS_WATCH(delta_acummulated)
                    }

                    returned_by_mapping_force = initial_tau_strength * calculation_area - k_unload * (delta_acummulated - delta_at_undamaged_peak);

                    if (returned_by_mapping_force < 0.0) {
                        returned_by_mapping_force = 0.0;
                    }

                    if (current_tangential_force_module) {
                        LocalElasticContactForce[0] = (returned_by_mapping_force / current_tangential_force_module) * LocalElasticContactForce[0];
                        LocalElasticContactForce[1] = (returned_by_mapping_force / current_tangential_force_module) * LocalElasticContactForce[1];
                    }

                    if (delta_acummulated) {
                        mKtUpdated = returned_by_mapping_force / delta_acummulated;
                    }

                    mDamage = 1.0 - returned_by_mapping_force/initial_tau_strength;
                    if (mDamage > damage_threshold_tolerance) {
                        failure_type = 2; // failure by shear
                    }
                    //
                    tau_strength = (1.0 - mDamage) * mTauZero;

                    if (contact_sigma >= 0) {
                        tau_strength += (1.0 - mDamage) * mInternalFriction * contact_sigma;
                    }
                    //
                }
            }
        } else {
            failure = 1;
            equiv_tg_of_fri_ang = 0.5 * (element1->GetTgOfFrictionAngle() + element2->GetTgOfFrictionAngle());
            maximum_frictional_shear_force = equiv_tg_of_fri_ang * LocalElasticContactForce[2];

            if (maximum_frictional_shear_force < 0.0) maximum_frictional_shear_force = 0.0;

            if (current_tangential_force_module > maximum_frictional_shear_force) {
                LocalElasticContactForce[0] = (maximum_frictional_shear_force / current_tangential_force_module) * LocalElasticContactForce[0];
                LocalElasticContactForce[1] = (maximum_frictional_shear_force / current_tangential_force_module) * LocalElasticContactForce[1];
                sliding = true;
                slide = 1;
            }
        }

        if (element1->Id() == 222222222222) {
            static std::ofstream tangential_forces_file("tangential_forces_damage.txt", std::ios_base::out | std::ios_base::app);
            tangential_forces_file << r_process_info[TIME] << " " << failure << " " << LocalElasticContactForce[0] << " " << tau_strength << " "
                                   << delta_acummulated << " " << returned_by_mapping_force << " " << mKtUpdated << " " << damage << " "
                                   << slide << " " << contact_sigma << " " << mDamage << " " << LocalDeltDisp[0] << " "
                                   << current_tangential_force_module << '\n';
            tangential_forces_file.flush();
        }

        KRATOS_CATCH("")
    }

} // namespace Kratos
