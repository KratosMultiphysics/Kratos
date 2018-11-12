// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_CamClay_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_CamClay::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_CamClay(*this));
        return p_clone;
    }

    void DEM_KDEM_CamClay::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_CamClay to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_KDEM_CamClay::CheckFailure(const int i_neighbour_count, SphericContinuumParticle* element1, SphericContinuumParticle* element2){

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

        if (failure_type == 0) {

            Matrix average_stress_tensor = ZeroMatrix(3,3);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) average_stress_tensor(i,j) = 0.5 * ((*(element1->mSymmStressTensor))(i,j) + (*(element2->mSymmStressTensor))(i,j));
            }

            Vector principal_stresses(3);
            noalias(principal_stresses) = AuxiliaryFunctions::EigenValuesDirectMethod(average_stress_tensor);

            Properties& element1_props = element1->GetProperties();
            Properties& element2_props = element2->GetProperties();

            // Preconsolidation pressure
            const double p_c = 0.5 * (element1_props[DEM_PRECONSOLIDATION_PRESSURE] + element2_props[DEM_PRECONSOLIDATION_PRESSURE]);
            
            // p and q computation
            const double p = 0.333333333333333333333 * (principal_stresses[0] + principal_stresses[1] + principal_stresses[2]);
            
            const double q = sqrt(0.5 * ((principal_stresses[0] - principal_stresses[1]) * (principal_stresses[0] - principal_stresses[1])
                                       + (principal_stresses[1] - principal_stresses[2]) * (principal_stresses[1] - principal_stresses[2])
                                       + (principal_stresses[2] - principal_stresses[0]) * (principal_stresses[2] - principal_stresses[0])));

            // slope of the straight line function
            const double M = 0.5 * (element1_props[DEM_M_CAMCLAY_SLOPE] + element2_props[DEM_M_CAMCLAY_SLOPE]);;
            
            // straight line function value
            const double straight_line_function_value = M * p;

            // ellipsoid function value
            const double ellipsoid_function_value = M * M * p * (p - p_c) + q * q;
            
            // choosing the necessary branch
            const double cam_clay_function_value = straight_line_function_value < ellipsoid_function_value ? straight_line_function_value : ellipsoid_function_value;

            if (cam_clay_function_value > 0) failure_type = 4;
        }
    }

    void DEM_KDEM_CamClay::CalculateNormalForces(double LocalElasticContactForce[3],
            const double kn_el,
            double equiv_young,
            double indentation,
            double calculation_area,
            double& acumulated_damage,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            int time_steps) {

        KRATOS_TRY

        //Firstly, we check that the bond is not broken (it can break in any state of forces or indentations, because breakage depends on the stress tensor)
        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

        if (indentation >= 0.0) { //COMPRESSION This response is the same for broken or intact bonds!
            LocalElasticContactForce[2] = kn_el * indentation;
        } else {
            if (failure_type > 0) {
                LocalElasticContactForce[2] = 0.0;
            } else {
                LocalElasticContactForce[2] = kn_el * indentation;
            }
        }

        KRATOS_CATCH("")
    }

    void DEM_KDEM_CamClay::CalculateTangentialForces(double OldLocalElasticContactForce[3],
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

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];
        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - kt_el * LocalDeltDisp[0]; // 0: first tangential
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - kt_el * LocalDeltDisp[1]; // 1: second tangential

        if (failure_type == 0) {
            if (r_process_info[SHEAR_STRAIN_PARALLEL_TO_BOND_OPTION]) {
                AddContributionOfShearStrainParallelToBond(OldLocalElasticContactForce,
                                                           LocalElasticExtraContactForce,
                                                           element1->mNeighbourElasticExtraContactForces[i_neighbour_count],
                                                           LocalCoordSystem,
                                                           kt_el, calculation_area,  element1, element2);
            }
        } else {
            LocalElasticExtraContactForce[0] = 0.0;
            LocalElasticExtraContactForce[1] = 0.0;

            double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
            const double equiv_tg_of_fri_ang = 0.5 * (element1->GetTgOfFrictionAngle() + element2->GetTgOfFrictionAngle());
            double Frictional_ShearForceMax = equiv_tg_of_fri_ang * LocalElasticContactForce[2];

            if (Frictional_ShearForceMax < 0.0) Frictional_ShearForceMax = 0.0;

            if ((ShearForceNow > Frictional_ShearForceMax) && (ShearForceNow != 0.0)) {
                LocalElasticContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[0];
                LocalElasticContactForce[1] = (Frictional_ShearForceMax / ShearForceNow) * LocalElasticContactForce[1];
                sliding = true;
            }
        }

        KRATOS_CATCH("")
    }
    
    double DEM_KDEM_CamClay::LocalMaxSearchDistance(const int i, SphericContinuumParticle* element1, SphericContinuumParticle* element2) {

        Properties& element1_props = element1->GetProperties();
        Properties& element2_props = element2->GetProperties();
        const double mean_preconsolidation_pressure = 1e6 * 0.5*(element1_props[DEM_PRECONSOLIDATION_PRESSURE] + element2_props[DEM_PRECONSOLIDATION_PRESSURE]);
        
        // calculation of equivalent young modulus
        const double myYoung = element1->GetYoung();
        const double other_young = element2->GetYoung();
        const double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);

        const double my_radius = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        
        double calculation_area = 0;
        Vector& vector_of_contact_areas = element1->GetValue(NEIGHBOURS_CONTACT_AREAS);
        GetContactArea(my_radius, other_radius, vector_of_contact_areas, i, calculation_area);
        
        const double radius_sum = my_radius + other_radius;
        const double initial_delta = element1->GetInitialDelta(i);
        const double initial_dist = radius_sum - initial_delta;

        const double kn_el = equiv_young * calculation_area / initial_dist;
        const double max_normal_force = mean_preconsolidation_pressure * calculation_area;

        const double max_local_distance_by_force = max_normal_force / kn_el;
        const double max_local_distance_by_radius = 0.05 * radius_sum;
        
        // avoid error in special cases with too high tensile
        const double max_local_distance = max_local_distance_by_force > max_local_distance_by_radius? max_local_distance_by_radius : max_local_distance_by_force;
                        
        return max_local_distance;
    }

    bool DEM_KDEM_CamClay::CheckRequirementsOfStressTensor() { return true;}

} // namespace Kratos
