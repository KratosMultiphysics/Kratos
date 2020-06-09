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

    void DEM_KDEM_CamClay::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_CamClay to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_KDEM_CamClay::CheckFailure(const int i_neighbour_count, SphericContinuumParticle* element1, SphericContinuumParticle* element2){

        int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

        if (failure_type == 0) {

            BoundedMatrix<double, 3, 3> average_stress_tensor = ZeroMatrix(3,3);
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

    double DEM_KDEM_CamClay::LocalMaxSearchDistance(const int i, SphericContinuumParticle* element1, SphericContinuumParticle* element2) {

        //Properties& element1_props = element1->GetProperties();
        //Properties& element2_props = element2->GetProperties();
        //const double mean_preconsolidation_pressure = 0.5*(element1_props[DEM_PRECONSOLIDATION_PRESSURE] + element2_props[DEM_PRECONSOLIDATION_PRESSURE]);

        BoundedMatrix<double, 3, 3> average_stress_tensor = ZeroMatrix(3,3);
        for (unsigned i = 0; i < 3; i++) {
            for (unsigned j = 0; j < 3; j++) {
                average_stress_tensor(i,j) = 0.5 * ((*(element1->mSymmStressTensor))(i,j) + (*(element2->mSymmStressTensor))(i,j));
            }
        }

        Vector principal_stresses(3);
        noalias(principal_stresses) = AuxiliaryFunctions::EigenValuesDirectMethod(average_stress_tensor);
        const double max_stress = *std::max_element(principal_stresses.begin(), principal_stresses.end());

        const double myYoung = element1->GetYoung();
        const double other_young = element2->GetYoung();
        const double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
        const double my_radius = element1->GetRadius();
        const double other_radius = element2->GetRadius();

        double calculation_area = 0.0;
        Vector& vector_of_contact_areas = element1->GetValue(NEIGHBOURS_CONTACT_AREAS);
        GetContactArea(my_radius, other_radius, vector_of_contact_areas, i, calculation_area);

        const double radius_sum = my_radius + other_radius;
        const double initial_delta = element1->GetInitialDelta(i);
        const double initial_dist = radius_sum - initial_delta;
        const double kn_el = equiv_young * calculation_area / initial_dist;
        //const double max_normal_force = mean_preconsolidation_pressure * calculation_area;
        const double max_normal_force = max_stress * calculation_area;
        const double max_local_distance_by_force = max_normal_force / kn_el;
        const double max_local_distance_by_radius = 0.05 * radius_sum;

        // avoid error in special cases with too high tensile
        const double max_local_distance = max_local_distance_by_force > max_local_distance_by_radius? max_local_distance_by_radius : max_local_distance_by_force;

        return max_local_distance;
    }
} // namespace Kratos
