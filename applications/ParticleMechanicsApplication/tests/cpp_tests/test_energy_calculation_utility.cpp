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
#include "particle_mechanics_application_variables.h"
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
        auto pElement = rModelPart.CreateNewElement("UpdatedLagrangian3D4N", 1, {{1, 2, 3, 4}}, p_elem_prop);

        // For potential energy
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 0.0;
        mp_coordinate[1] = 0.5;
        mp_coordinate[2] = 0.0;
        array_1d<double, 3> volume_acceleration;
        volume_acceleration[0] = 0.0;
        volume_acceleration[1] = -9.8;
        volume_acceleration[2] = 0.0;
        pElement->SetValue(MP_COORD, mp_coordinate);
        pElement->SetValue(MP_MASS, 1.5);
        pElement->SetValue(MP_VOLUME_ACCELERATION, volume_acceleration);

        // For kinetic energy
        array_1d<double, 3> velocity;
        velocity[0] = 1.0;
        velocity[1] = 2.0;
        velocity[2] = 3.0;
        pElement->SetValue(MP_VELOCITY, velocity);

        // For strain energy
        array_1d<double, 6> mp_cauchy_stress;
        mp_cauchy_stress[0] = 1.0;
        mp_cauchy_stress[1] = 2.0;
        mp_cauchy_stress[2] = 3.0;
        mp_cauchy_stress[3] = 4.0;
        mp_cauchy_stress[4] = 5.0;
        mp_cauchy_stress[5] = 6.0;
        array_1d<double, 6> mp_strain;
        mp_strain[0] = 0.1;
        mp_strain[1] = 0.2;
        mp_strain[2] = 0.3;
        mp_strain[3] = 0.4;
        mp_strain[4] = 0.5;
        mp_strain[5] = 0.6;
        pElement->SetValue(MP_VOLUME, 2.5);
        pElement->SetValue(MP_CAUCHY_STRESS_VECTOR, mp_cauchy_stress);
        pElement->SetValue(MP_ALMANSI_STRAIN_VECTOR, mp_strain);

    }

    /**
    * Check whether the calculation of energy are okay
    */
    KRATOS_TEST_CASE_IN_SUITE(ParticleTotalEnergyCalculation, KratosParticleMechanicsFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        PrepareModelPart(r_model_part);

        // Check energy
        MPMEnergyCalculationUtility::CalculateTotalEnergy(r_model_part);

        const std::size_t element_id = 1;
        const double & r_MP_PotentialEnergy = r_model_part.pGetElement(element_id)->GetValue(MP_POTENTIAL_ENERGY);
        const double & r_MP_KineticEnergy   = r_model_part.pGetElement(element_id)->GetValue(MP_KINETIC_ENERGY);
        const double & r_MP_StrainEnergy    = r_model_part.pGetElement(element_id)->GetValue(MP_STRAIN_ENERGY);
        const double & r_MP_TotalEnergy     = r_model_part.pGetElement(element_id)->GetValue(MP_TOTAL_ENERGY);
        KRATOS_CHECK_NEAR(r_MP_PotentialEnergy, 7.35 , 1e-6);
        KRATOS_CHECK_NEAR(r_MP_KineticEnergy  ,10.50 , 1e-6);
        KRATOS_CHECK_NEAR(r_MP_StrainEnergy   ,11.375, 1e-6);
        KRATOS_CHECK_NEAR(r_MP_TotalEnergy    ,29.225, 1e-6);
    }

} // namespace Testing
} // namespace Kratos
