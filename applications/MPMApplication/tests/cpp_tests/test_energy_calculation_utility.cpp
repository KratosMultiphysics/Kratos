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
#include "containers/model.h"
#include "custom_utilities/mpm_energy_calculation_utility.h"
#include "mpm_application_variables.h"

namespace Kratos::Testing
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
        auto pElement1 = rModelPart.CreateNewElement("MPMUpdatedLagrangian3D4N", 1, {{1, 2, 3, 4}}, p_elem_prop);
        auto pElement2 = rModelPart.CreateNewElement("MPMUpdatedLagrangian3D4N", 2, {{1, 2, 3, 4}}, p_elem_prop);
        auto pElement3 = rModelPart.CreateNewElement("MPMUpdatedLagrangian3D4N", 3, {{1, 2, 3, 4}}, p_elem_prop);

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // coordinates
        array_1d<double, 3> mp_coordinate1{0.0, 0.5, 0.0};
        pElement1->SetValuesOnIntegrationPoints(MP_COORD, { mp_coordinate1 }, r_current_process_info);
        array_1d<double, 3> mp_coordinate2{0.3, 1.5, 0.0};
        pElement2->SetValuesOnIntegrationPoints(MP_COORD, { mp_coordinate2 }, r_current_process_info);
        array_1d<double, 3> mp_coordinate3{1.1, 0.3, 1.2};
        pElement3->SetValuesOnIntegrationPoints(MP_COORD, { mp_coordinate3 }, r_current_process_info);

        // volume acceleration
        array_1d<double, 3> volume_acceleration1{0.0, -9.8, 0.0};
        pElement1->SetValuesOnIntegrationPoints(MP_VOLUME_ACCELERATION, { volume_acceleration1 }, r_current_process_info);
        array_1d<double, 3> volume_acceleration2{-9.8, 0.0, 0.0};
        pElement2->SetValuesOnIntegrationPoints(MP_VOLUME_ACCELERATION, { volume_acceleration2 }, r_current_process_info);
        array_1d<double, 3> volume_acceleration3{0.0, 2.0, 1.0};
        pElement3->SetValuesOnIntegrationPoints(MP_VOLUME_ACCELERATION, { volume_acceleration3 }, r_current_process_info);

        // mass
        double mp_mass1{ 1.5 };
        pElement1->SetValuesOnIntegrationPoints(MP_MASS, { mp_mass1 }, r_current_process_info);
        double mp_mass2{ 3.0 };
        pElement2->SetValuesOnIntegrationPoints(MP_MASS, { mp_mass2 }, r_current_process_info);
        double mp_mass3{ 2.2 };
        pElement3->SetValuesOnIntegrationPoints(MP_MASS, { mp_mass3 }, r_current_process_info);

        // For kinetic energy
        array_1d<double, 3> velocity1{1.0, 2.0, 3.0};
        pElement1->SetValuesOnIntegrationPoints(MP_VELOCITY, { velocity1 }, r_current_process_info);
        array_1d<double, 3> velocity2{1.3, 1.5, 1.0};
        pElement2->SetValuesOnIntegrationPoints(MP_VELOCITY, { velocity2 }, r_current_process_info);
        array_1d<double, 3> velocity3{0.9, 0.5, 0.6};
        pElement3->SetValuesOnIntegrationPoints(MP_VELOCITY, { velocity3 }, r_current_process_info);

        // For strain energy
        double mp_volume1{ 2.5 };
        pElement1->SetValuesOnIntegrationPoints(MP_VOLUME, { mp_volume1 }, r_current_process_info);
        double mp_volume2{ 0.5 };
        pElement2->SetValuesOnIntegrationPoints(MP_VOLUME, { mp_volume2 }, r_current_process_info);
        double mp_volume3{ 1.8 };
        pElement3->SetValuesOnIntegrationPoints(MP_VOLUME, { mp_volume3 }, r_current_process_info);

        Vector mp_cauchy_stress1 = ZeroVector(6);
        mp_cauchy_stress1[0] = 1.0;
        mp_cauchy_stress1[1] = 2.0;
        mp_cauchy_stress1[2] = 3.0;
        mp_cauchy_stress1[3] = 4.0;
        mp_cauchy_stress1[4] = 5.0;
        mp_cauchy_stress1[5] = 6.0;
        pElement1->SetValuesOnIntegrationPoints(MP_CAUCHY_STRESS_VECTOR, { mp_cauchy_stress1 }, r_current_process_info);
        Vector mp_cauchy_stress2 = ZeroVector(6);
        mp_cauchy_stress2[0] = 0.2;
        mp_cauchy_stress2[1] = 0.1;
        mp_cauchy_stress2[2] = -0.8;
        mp_cauchy_stress2[3] = -0.3;
        mp_cauchy_stress2[4] = 0.2;
        mp_cauchy_stress2[5] = 1.2;
        pElement2->SetValuesOnIntegrationPoints(MP_CAUCHY_STRESS_VECTOR, { mp_cauchy_stress2 }, r_current_process_info);
        Vector mp_cauchy_stress3 = ZeroVector(6);
        mp_cauchy_stress3[0] = -3.0;
        mp_cauchy_stress3[1] = 0.6;
        mp_cauchy_stress3[2] = -0.1;
        mp_cauchy_stress3[3] = 0.0;
        mp_cauchy_stress3[4] = 0.1;
        mp_cauchy_stress3[5] = 0.8;
        pElement3->SetValuesOnIntegrationPoints(MP_CAUCHY_STRESS_VECTOR, { mp_cauchy_stress3 }, r_current_process_info);

        Vector mp_almansi_strain1 = ZeroVector(6);
        mp_almansi_strain1[0] = 0.1;
        mp_almansi_strain1[1] = 0.2;
        mp_almansi_strain1[2] = 0.3;
        mp_almansi_strain1[3] = 0.4;
        mp_almansi_strain1[4] = 0.5;
        mp_almansi_strain1[5] = 0.6;
        pElement1->SetValuesOnIntegrationPoints(MP_ALMANSI_STRAIN_VECTOR, { mp_almansi_strain1 }, r_current_process_info);
        Vector mp_almansi_strain2 = ZeroVector(6);
        mp_almansi_strain2[0] = -0.3;
        mp_almansi_strain2[1] = -0.7;
        mp_almansi_strain2[2] = 0.2;
        mp_almansi_strain2[3] = 1.9;
        mp_almansi_strain2[4] = 0.1;
        mp_almansi_strain2[5] = 1.6;
        pElement2->SetValuesOnIntegrationPoints(MP_ALMANSI_STRAIN_VECTOR, { mp_almansi_strain2 }, r_current_process_info);
        Vector mp_almansi_strain3 = ZeroVector(6);
        mp_almansi_strain3[0] = 0.0;
        mp_almansi_strain3[1] = 0.5;
        mp_almansi_strain3[2] = -0.3;
        mp_almansi_strain3[3] = 2.4;
        mp_almansi_strain3[4] = -0.6;
        mp_almansi_strain3[5] = -1.6;
        pElement3->SetValuesOnIntegrationPoints(MP_ALMANSI_STRAIN_VECTOR, { mp_almansi_strain3 }, r_current_process_info);
    }

    /**
    * Check whether the calculation of energy are okay
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMTotalEnergyCalculationElement, KratosMPMFastSuite)
    {
        KRATOS_WATCH("")
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        PrepareModelPart(r_model_part);

        const ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();

        std::vector<double> r_MP_PotentialEnergy(1);
        std::vector<double> r_MP_KineticEnergy(1);
        std::vector<double> r_MP_StrainEnergy(1);
        std::vector<double> r_MP_TotalEnergy(1);

        r_model_part.pGetElement(1)->CalculateOnIntegrationPoints(MP_POTENTIAL_ENERGY, r_MP_PotentialEnergy, r_current_process_info);
        r_model_part.pGetElement(1)->CalculateOnIntegrationPoints(MP_KINETIC_ENERGY, r_MP_KineticEnergy, r_current_process_info);
        r_model_part.pGetElement(1)->CalculateOnIntegrationPoints(MP_STRAIN_ENERGY, r_MP_StrainEnergy, r_current_process_info);
        r_model_part.pGetElement(1)->CalculateOnIntegrationPoints(MP_TOTAL_ENERGY, r_MP_TotalEnergy, r_current_process_info);

        KRATOS_EXPECT_NEAR(r_MP_PotentialEnergy[0], 7.35  , 1e-6);
        KRATOS_EXPECT_NEAR(r_MP_KineticEnergy[0]  , 10.50 , 1e-6);
        KRATOS_EXPECT_NEAR(r_MP_StrainEnergy[0]   , 11.375, 1e-6);
        KRATOS_EXPECT_NEAR(r_MP_TotalEnergy[0]    , 29.225, 1e-6);

        r_model_part.pGetElement(2)->CalculateOnIntegrationPoints(MP_POTENTIAL_ENERGY, r_MP_PotentialEnergy, r_current_process_info);
        r_model_part.pGetElement(2)->CalculateOnIntegrationPoints(MP_KINETIC_ENERGY, r_MP_KineticEnergy, r_current_process_info);
        r_model_part.pGetElement(2)->CalculateOnIntegrationPoints(MP_STRAIN_ENERGY, r_MP_StrainEnergy, r_current_process_info);
        r_model_part.pGetElement(2)->CalculateOnIntegrationPoints(MP_TOTAL_ENERGY, r_MP_TotalEnergy, r_current_process_info);

        KRATOS_EXPECT_NEAR(r_MP_PotentialEnergy[0],  8.82 , 1e-6);
        KRATOS_EXPECT_NEAR(r_MP_KineticEnergy[0]  ,  7.41 , 1e-6);
        KRATOS_EXPECT_NEAR(r_MP_StrainEnergy[0]   ,  0.27 , 1e-6);
        KRATOS_EXPECT_NEAR(r_MP_TotalEnergy[0]    , 16.50 , 1e-6);

        r_model_part.pGetElement(3)->CalculateOnIntegrationPoints(MP_POTENTIAL_ENERGY, r_MP_PotentialEnergy, r_current_process_info);
        r_model_part.pGetElement(3)->CalculateOnIntegrationPoints(MP_KINETIC_ENERGY, r_MP_KineticEnergy, r_current_process_info);
        r_model_part.pGetElement(3)->CalculateOnIntegrationPoints(MP_STRAIN_ENERGY, r_MP_StrainEnergy, r_current_process_info);
        r_model_part.pGetElement(3)->CalculateOnIntegrationPoints(MP_TOTAL_ENERGY, r_MP_TotalEnergy, r_current_process_info);

        KRATOS_EXPECT_NEAR(r_MP_PotentialEnergy[0],  3.96  , 1e-6);
        KRATOS_EXPECT_NEAR(r_MP_KineticEnergy[0]  ,  1.562 , 1e-6);
        KRATOS_EXPECT_NEAR(r_MP_StrainEnergy[0]   , -0.909 , 1e-6);
        KRATOS_EXPECT_NEAR(r_MP_TotalEnergy[0]    ,  4.613 , 1e-6);

        // Check energy
        const double potential_energy_model_part = MPMEnergyCalculationUtility().CalculatePotentialEnergy(r_model_part);
        const double kinetic_energy_model_part = MPMEnergyCalculationUtility().CalculateKineticEnergy(r_model_part);
        const double strain_energy_model_part = MPMEnergyCalculationUtility().CalculateStrainEnergy(r_model_part);
        const double total_energy_model_part = MPMEnergyCalculationUtility().CalculateTotalEnergy(r_model_part);

        KRATOS_EXPECT_NEAR(potential_energy_model_part, 20.13  , 1e-6);
        KRATOS_EXPECT_NEAR(kinetic_energy_model_part,   19.472 , 1e-6);
        KRATOS_EXPECT_NEAR(strain_energy_model_part,    10.736 , 1e-6);
        KRATOS_EXPECT_NEAR(total_energy_model_part,     50.338 , 1e-6);

        double p_energy{0};
        double k_energy{0};
        double s_energy{0};
        double t_energy{0};

        MPMEnergyCalculationUtility().CalculateTotalEnergy(r_model_part, p_energy, k_energy, s_energy, t_energy);

        KRATOS_EXPECT_NEAR(p_energy, 20.13  , 1e-6);
        KRATOS_EXPECT_NEAR(k_energy, 19.472 , 1e-6);
        KRATOS_EXPECT_NEAR(s_energy, 10.736 , 1e-6);
        KRATOS_EXPECT_NEAR(t_energy, 50.338 , 1e-6);

    }

} // namespace Kratos::Testing
