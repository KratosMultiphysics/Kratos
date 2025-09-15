//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "custom_utilities/mpm_energy_calculation_utility.h"
#include "mpm_application_variables.h"
#include "containers/model.h"

namespace Kratos
{
namespace Testing
{
    void PrepareModelPart(ModelPart& rModelPart)
    {
        // Nodes
        auto p_node_1 = rModelPart.CreateNewNode( 1,  0.0 ,  0.0 , 0.0);
        auto p_node_2 = rModelPart.CreateNewNode( 2,  1.0 ,  0.0 , 0.0);
        auto p_node_3 = rModelPart.CreateNewNode( 3,  1.0 ,  1.0 , 0.0);
        auto p_node_4 = rModelPart.CreateNewNode( 4,  0.0 ,  1.0 , 0.0);

        // Properties
        Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

        // Elements
        auto pElement = rModelPart.CreateNewElement("MPMUpdatedLagrangian3D4N", 1, {{1, 2, 3, 4}}, p_elem_prop);

        // For potential energy
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 0.0;
        mp_coordinate[1] = 0.5;
        mp_coordinate[2] = 0.0;
        array_1d<double, 3> volume_acceleration;
        volume_acceleration[0] = 0.0;
        volume_acceleration[1] = -9.8;
        volume_acceleration[2] = 0.0;
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        pElement->SetValuesOnIntegrationPoints(MP_COORD, { mp_coordinate }, r_current_process_info);
        std::vector<double> mp_mass_vector = { 1.5 };
        pElement->SetValuesOnIntegrationPoints(MP_MASS, mp_mass_vector, r_current_process_info);
        pElement->SetValuesOnIntegrationPoints(MP_VOLUME_ACCELERATION, { volume_acceleration }, r_current_process_info);

        // For kinetic energy
        array_1d<double, 3> velocity;
        velocity[0] = 1.0;
        velocity[1] = 2.0;
        velocity[2] = 3.0;
        pElement->SetValuesOnIntegrationPoints(MP_VELOCITY, { velocity }, r_current_process_info);

        // For strain energy
        //Vector mp_cauchy_stress1 = ;
        Vector mp_cauchy_stress = ZeroVector(6);
        mp_cauchy_stress[0] = 1.0;
        mp_cauchy_stress[1] = 2.0;
        mp_cauchy_stress[2] = 3.0;
        mp_cauchy_stress[3] = 4.0;
        mp_cauchy_stress[4] = 5.0;
        mp_cauchy_stress[5] = 6.0;
        Vector mp_strain = ZeroVector(6);
        mp_strain[0] = 0.1;
        mp_strain[1] = 0.2;
        mp_strain[2] = 0.3;
        mp_strain[3] = 0.4;
        mp_strain[4] = 0.5;
        mp_strain[5] = 0.6;
        std::vector<double> mp_volume_vector = { 2.5 };
        pElement->SetValuesOnIntegrationPoints(MP_VOLUME, mp_volume_vector, r_current_process_info);
        std::vector<Vector> mp_cauchy_stress_vector = { mp_cauchy_stress };
        pElement->SetValuesOnIntegrationPoints(MP_CAUCHY_STRESS_VECTOR, mp_cauchy_stress_vector, r_current_process_info);
        std::vector<Vector> mp_almansi_strain_vector = { mp_strain };
        pElement->SetValuesOnIntegrationPoints(MP_ALMANSI_STRAIN_VECTOR, mp_almansi_strain_vector, r_current_process_info);
    }

    /**
    * Check whether the calculation of energy are okay
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMTotalEnergyCalculation, KratosMPMFastSuite)
    {
        KRATOS_WATCH("")
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        PrepareModelPart(r_model_part);

        // Check energy
        MPMEnergyCalculationUtility::CalculateTotalEnergy(r_model_part);

        const std::size_t element_id = 1;

        std::vector<double> r_MP_PotentialEnergy(1);
        const ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();

        r_model_part.pGetElement(element_id)->CalculateOnIntegrationPoints(MP_POTENTIAL_ENERGY, r_MP_PotentialEnergy, r_current_process_info);
        
        std::vector<double> r_MP_KineticEnergy(1);
        r_model_part.pGetElement(element_id)->CalculateOnIntegrationPoints(MP_KINETIC_ENERGY, r_MP_KineticEnergy, r_current_process_info);

        std::vector<double> r_MP_StrainEnergy(1);
        r_model_part.pGetElement(element_id)->CalculateOnIntegrationPoints(MP_STRAIN_ENERGY, r_MP_StrainEnergy, r_current_process_info);

        std::vector<double> r_MP_TotalEnergy(1);
        r_model_part.pGetElement(element_id)->CalculateOnIntegrationPoints(MP_TOTAL_ENERGY, r_MP_TotalEnergy, r_current_process_info);

        KRATOS_EXPECT_NEAR(r_MP_PotentialEnergy[0], 7.35  , 1e-6);
        KRATOS_EXPECT_NEAR(r_MP_KineticEnergy[0]  , 10.50 , 1e-6);
        KRATOS_EXPECT_NEAR(r_MP_StrainEnergy[0]   , 11.375, 1e-6);
        KRATOS_EXPECT_NEAR(r_MP_TotalEnergy[0]    , 29.225, 1e-6);
    }

} // namespace Testing
} // namespace Kratos
