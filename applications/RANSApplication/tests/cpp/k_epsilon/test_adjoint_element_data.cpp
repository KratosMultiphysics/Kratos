//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/cfd_variables.h"
#include "testing/testing.h"

// Application includes
#include "../stabilization_method_test_utilities.h"
#include "custom_elements/convection_diffusion_reaction_element_data/k_epsilon/epsilon_adjoint_element_data.h"
#include "custom_elements/convection_diffusion_reaction_element_data/k_epsilon/epsilon_element_data.h"
#include "custom_elements/convection_diffusion_reaction_element_data/k_epsilon/k_adjoint_element_data.h"
#include "custom_elements/convection_diffusion_reaction_element_data/k_epsilon/k_element_data.h"
#include "custom_processes/rans_nut_k_epsilon_update_process.h"
#include "custom_utilities/adjoint_test_utilities.h"
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"
#include "test_utilities.h"

namespace Kratos
{
namespace Testing
{
namespace
{
void KEpsilon_AddVariables(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_3);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_3);
    rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_1);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_3);
    rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_3);
}

void KEpsilon_SetVariables(ModelPart& rModelPart)
{
    using namespace RansApplicationTestUtilities;

    RandomFillNodalHistoricalVariable(rModelPart, VELOCITY, -10.0, 10.0);
    RandomFillNodalHistoricalVariable(rModelPart, KINEMATIC_VISCOSITY, 1e-5, 1e-3);
    RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    RandomFillNodalHistoricalVariable(
        rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE, 1.0, 100.0);
    RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 1.0, 50.0);
    RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE_2, 1.0, 50.0);
    RandomFillNodalHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_1, 1.0, 10.0);
    RandomFillNodalHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_2, 1.0, 10.0);

    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 0.5);
    r_process_info.SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.5);
    r_process_info.SetValue(TURBULENCE_RANS_C1, 2.5);
    r_process_info.SetValue(TURBULENCE_RANS_C2, 3.5);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 2.1);
}

void KEpsilon_UpdateVariables(ModelPart& rModelPart)
{
    RansNutKEpsilonUpdateProcess nut_update(rModelPart.GetModel(),
                                            rModelPart.Name(), 2.1, 1e-12, 0);
    nut_update.Execute();
}

} // namespace

/// KAdjointElementData CalculateEffectiveKinematicViscosityDerivatives tests
KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateEffectiveKinematicViscosityDerivatives_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_KINETIC_ENERGY, &primal_type::CalculateEffectiveKinematicViscosity,
        &adjoint_type::CalculateEffectiveKinematicViscosityDerivatives, 1, 1e-6,
        1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateEffectiveKinematicViscosityDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, TURBULENT_ENERGY_DISSIPATION_RATE,
        &primal_type::CalculateEffectiveKinematicViscosity,
        &adjoint_type::CalculateEffectiveKinematicViscosityDerivatives, 1, 1e-6,
        1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateEffectiveKinematicViscosityDerivatives_VELOCITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        VELOCITY, &primal_type::CalculateEffectiveKinematicViscosity,
        &adjoint_type::CalculateEffectiveKinematicViscosityDerivatives, 1, 1e-6,
        1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateEffectiveKinematicViscosityDerivatives_SHAPE_SENSITIVITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, &primal_type::CalculateEffectiveKinematicViscosity,
        &adjoint_type::CalculateEffectiveKinematicViscosityShapeDerivatives, 1,
        1e-6, 1e-6, 1e-9);
}

/// KAdjointElementData CalculateReactionTermDerivatives tests
KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateReactionTermDerivatives_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_KINETIC_ENERGY, &primal_type::CalculateReactionTerm,
        &adjoint_type::CalculateReactionTermDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateReactionTermDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_ENERGY_DISSIPATION_RATE, &primal_type::CalculateReactionTerm,
        &adjoint_type::CalculateReactionTermDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateReactionTermDerivatives_VELOCITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, VELOCITY, &primal_type::CalculateReactionTerm,
        &adjoint_type::CalculateReactionTermDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateReactionTermDerivatives_SHAPE_SENSITIVITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, &primal_type::CalculateReactionTerm,
        &adjoint_type::CalculateReactionTermShapeDerivatives, 1, 1e-7, 1e-6, 1e-9);
}

/// KAdjointElementData CalculateSourceTermDerivatives tests
KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateSourceTermDerivatives_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_KINETIC_ENERGY, &primal_type::CalculateSourceTerm,
        &adjoint_type::CalculateSourceTermDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateSourceTermDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_ENERGY_DISSIPATION_RATE, &primal_type::CalculateSourceTerm,
        &adjoint_type::CalculateSourceTermDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateSourceTermDerivatives_VELOCITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, VELOCITY, &primal_type::CalculateSourceTerm,
        &adjoint_type::CalculateSourceTermDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateSourceTermDerivatives_SHAPE_SENSITIVITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, &primal_type::CalculateSourceTerm,
        &adjoint_type::CalculateSourceTermShapeDerivatives, 1, 1e-7, 1e-6, 1e-9);
}

/// KAdjointElementData CalculateEffectiveVelocityDerivatives tests
KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateEffectiveVelocityDerivatives_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_KINETIC_ENERGY, &primal_type::CalculateEffectiveVelocity,
        &adjoint_type::CalculateEffectiveVelocityDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateEffectiveVelocityDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_ENERGY_DISSIPATION_RATE, &primal_type::CalculateEffectiveVelocity,
        &adjoint_type::CalculateEffectiveVelocityDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateEffectiveVelocityDerivatives_VELOCITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, VELOCITY, &primal_type::CalculateEffectiveVelocity,
        &adjoint_type::CalculateEffectiveVelocityDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_KAdjointElementData_CalculateEffectiveVelocityDerivatives_SHAPE_SENSITIVITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::KElementData<2>;
    using adjoint_type = KEpsilonElementData::KAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, &primal_type::CalculateEffectiveVelocity,
        &adjoint_type::CalculateEffectiveVelocityShapeDerivatives, 1, 1e-7, 1e-6, 1e-9);
}

/// EpsilonAdjointElementData CalculateEffectiveKinematicViscosityDerivatives tests
KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateEffectiveKinematicViscosityDerivatives_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_KINETIC_ENERGY, &primal_type::CalculateEffectiveKinematicViscosity,
        &adjoint_type::CalculateEffectiveKinematicViscosityDerivatives, 1, 1e-6,
        1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateEffectiveKinematicViscosityDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, TURBULENT_ENERGY_DISSIPATION_RATE,
        &primal_type::CalculateEffectiveKinematicViscosity,
        &adjoint_type::CalculateEffectiveKinematicViscosityDerivatives, 1, 1e-6,
        1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateEffectiveKinematicViscosityDerivatives_VELOCITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        VELOCITY, &primal_type::CalculateEffectiveKinematicViscosity,
        &adjoint_type::CalculateEffectiveKinematicViscosityDerivatives, 1, 1e-6,
        1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateEffectiveKinematicViscosityDerivatives_SHAPE_SENSITIVITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, &primal_type::CalculateEffectiveKinematicViscosity,
        &adjoint_type::CalculateEffectiveKinematicViscosityShapeDerivatives, 1,
        1e-6, 1e-6, 1e-9);
}

/// EpsilonAdjointElementData CalculateReactionTermDerivatives tests
KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateReactionTermDerivatives_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_KINETIC_ENERGY, &primal_type::CalculateReactionTerm,
        &adjoint_type::CalculateReactionTermDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateReactionTermDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_ENERGY_DISSIPATION_RATE, &primal_type::CalculateReactionTerm,
        &adjoint_type::CalculateReactionTermDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateReactionTermDerivatives_VELOCITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, VELOCITY, &primal_type::CalculateReactionTerm,
        &adjoint_type::CalculateReactionTermDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateReactionTermDerivatives_SHAPE_SENSITIVITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, &primal_type::CalculateReactionTerm,
        &adjoint_type::CalculateReactionTermShapeDerivatives, 1, 1e-7, 1e-6, 1e-9);
}

/// EpsilonAdjointElementData CalculateSourceTermDerivatives tests
KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateSourceTermDerivatives_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_KINETIC_ENERGY, &primal_type::CalculateSourceTerm,
        &adjoint_type::CalculateSourceTermDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateSourceTermDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_ENERGY_DISSIPATION_RATE, &primal_type::CalculateSourceTerm,
        &adjoint_type::CalculateSourceTermDerivatives, 1, 1e-3, 1e-6, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateSourceTermDerivatives_VELOCITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, VELOCITY, &primal_type::CalculateSourceTerm,
        &adjoint_type::CalculateSourceTermDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateSourceTermDerivatives_SHAPE_SENSITIVITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, &primal_type::CalculateSourceTerm,
        &adjoint_type::CalculateSourceTermShapeDerivatives, 1, 1e-7, 1e-6, 1e-9);
}

/// EpsilonAdjointElementData CalculateEffectiveVelocityDerivatives tests
KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateEffectiveVelocityDerivatives_TURBULENT_KINETIC_ENERGY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_KINETIC_ENERGY, &primal_type::CalculateEffectiveVelocity,
        &adjoint_type::CalculateEffectiveVelocityDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateEffectiveVelocityDerivatives_TURBULENT_ENERGY_DISSIPATION_RATE,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables, KEpsilon_UpdateVariables,
        TURBULENT_ENERGY_DISSIPATION_RATE, &primal_type::CalculateEffectiveVelocity,
        &adjoint_type::CalculateEffectiveVelocityDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateEffectiveVelocityDerivatives_VELOCITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, VELOCITY, &primal_type::CalculateEffectiveVelocity,
        &adjoint_type::CalculateEffectiveVelocityDerivatives, 1, 1e-6, 1e-6, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilon_EpsilonAdjointElementData_CalculateEffectiveVelocityDerivatives_SHAPE_SENSITIVITY,
                          KratosRansFastSuite)
{
    Model test_model;

    using primal_type = KEpsilonElementData::EpsilonElementData<2>;
    using adjoint_type = KEpsilonElementData::EpsilonAdjointElementData<2, 3>;

    RansApplicationTestUtilities::RunAdjointElementDataTest<primal_type, adjoint_type>(
        test_model, KEpsilon_AddVariables, KEpsilon_SetVariables,
        KEpsilon_UpdateVariables, &primal_type::CalculateEffectiveVelocity,
        &adjoint_type::CalculateEffectiveVelocityShapeDerivatives, 1, 1e-7, 1e-6, 1e-9);
}

} // namespace Testing
} // namespace Kratos.
