//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "../stabilization_method_test_utilities.h"
#include "custom_elements/data_containers/k_epsilon/k_adjoint_element_data.h"
#include "custom_elements/data_containers/k_epsilon/k_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_adjoint_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"
#include "includes/cfd_variables.h"
#include "test_utilities.h"
#include "custom_processes/rans_nut_k_epsilon_update_process.h"
#include "custom_utilities/rans_calculation_utilities.h"

#include "custom_utilities/adjoint_test_utilities.h"

namespace Kratos
{
namespace Testing
{

namespace
{
void KEpsilonAddVariables(ModelPart& rModelPart)
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

void KEpsilonSetVariables(ModelPart& rModelPart)
{
    using namespace RansApplicationTestUtilities;

    RandomFillNodalHistoricalVariable(rModelPart, VELOCITY, -10.0, 10.0);
    RandomFillNodalHistoricalVariable(rModelPart, KINEMATIC_VISCOSITY, 1e-5, 1e-3);
    RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE, 1.0, 100.0);
    RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 1.0, 50.0);
    RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE_2, 1.0, 50.0);
    RandomFillNodalHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_1, 1.0, 10.0);
    RandomFillNodalHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_2, 1.0, 10.0);

    auto& r_process_info = rModelPart.GetProcessInfo();
    r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 0.5);
    r_process_info.SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.5);
    r_process_info.SetValue(TURBULENCE_RANS_C1, 2.5);
    r_process_info.SetValue(TURBULENCE_RANS_C2, 3.5);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 2.1);
}

void KEpsilonUpdateVariables(ModelPart& rModelPart)
{
    RansNutKEpsilonUpdateProcess nut_update(rModelPart.GetModel(),
                                            rModelPart.Name(), 2.1, 1e-12, 0);
    nut_update.ExecuteAfterCouplingSolveStep();
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(RANSAdjointKEpsilonKEquationVelocityDerivatives, KratosRansFastSuite)
{
    Model model;

    RansApplicationTestUtilities::RunAdjointElementDataTest<
        KEpsilonElementData::KElementData<2>,
        KEpsilonAdjointElementData::KAdjointStateDerivatives<2, 3>::Data,
        KEpsilonAdjointElementData::KAdjointStateDerivatives<2, 3>::VelocityDerivatives>(
        model, KEpsilonAddVariables, KEpsilonSetVariables, KEpsilonUpdateVariables, 2, 1e-6, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RANSAdjointKEpsilonKEquationKDerivatives, KratosRansFastSuite)
{
    Model model;

    RansApplicationTestUtilities::RunAdjointElementDataTest<
        KEpsilonElementData::KElementData<2>,
        KEpsilonAdjointElementData::KAdjointStateDerivatives<2, 3>::Data,
        KEpsilonAdjointElementData::KAdjointStateDerivatives<2, 3>::KDerivatives>(
        model, KEpsilonAddVariables, KEpsilonSetVariables, KEpsilonUpdateVariables, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RANSAdjointKEpsilonKEquationEpsilonDerivatives, KratosRansFastSuite)
{
    Model model;

    RansApplicationTestUtilities::RunAdjointElementDataTest<
        KEpsilonElementData::KElementData<2>,
        KEpsilonAdjointElementData::KAdjointStateDerivatives<2, 3>::Data,
        KEpsilonAdjointElementData::KAdjointStateDerivatives<2, 3>::EpsilonDerivatives>(
        model, KEpsilonAddVariables, KEpsilonSetVariables, KEpsilonUpdateVariables, 2, 1e-6, 1e-4);
}

KRATOS_TEST_CASE_IN_SUITE(RANSAdjointKEpsilonKEquationShapeDerivatives, KratosRansFastSuite)
{
    Model model;

    RansApplicationTestUtilities::RunAdjointSensitivityDataTest<
        KEpsilonElementData::KElementData<2>,
        KEpsilonAdjointElementData::KAdjointShapeDerivatives<2, 3>::Data,
        KEpsilonAdjointElementData::KAdjointShapeDerivatives<2, 3>::ShapeDerivatives>(
        model, KEpsilonAddVariables, KEpsilonSetVariables, KEpsilonUpdateVariables, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RANSAdjointKEpsilonEpsilonEquationVelocityDerivatives, KratosRansFastSuite)
{
    Model model;

    RansApplicationTestUtilities::RunAdjointElementDataTest<
        KEpsilonElementData::EpsilonElementData<2>,
        KEpsilonAdjointElementData::EpsilonAdjointStateDerivatives<2, 3>::Data,
        KEpsilonAdjointElementData::EpsilonAdjointStateDerivatives<2, 3>::VelocityDerivatives>(
        model, KEpsilonAddVariables, KEpsilonSetVariables, KEpsilonUpdateVariables, 2, 1e-6, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RANSAdjointKEpsilonEpsilonEquationKDerivatives, KratosRansFastSuite)
{
    Model model;

    RansApplicationTestUtilities::RunAdjointElementDataTest<
        KEpsilonElementData::EpsilonElementData<2>,
        KEpsilonAdjointElementData::EpsilonAdjointStateDerivatives<2, 3>::Data,
        KEpsilonAdjointElementData::EpsilonAdjointStateDerivatives<2, 3>::KDerivatives>(
        model, KEpsilonAddVariables, KEpsilonSetVariables, KEpsilonUpdateVariables, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RANSAdjointKEpsilonEpsilonEquationEpsilonDerivatives, KratosRansFastSuite)
{
    Model model;

    RansApplicationTestUtilities::RunAdjointElementDataTest<
        KEpsilonElementData::EpsilonElementData<2>,
        KEpsilonAdjointElementData::EpsilonAdjointStateDerivatives<2, 3>::Data,
        KEpsilonAdjointElementData::EpsilonAdjointStateDerivatives<2, 3>::EpsilonDerivatives>(
        model, KEpsilonAddVariables, KEpsilonSetVariables, KEpsilonUpdateVariables, 2, 1e-6, 1e-4);
}

KRATOS_TEST_CASE_IN_SUITE(RANSAdjointKEpsilonEpsilonEquationShapeDerivatives, KratosRansFastSuite)
{
    Model model;

    RansApplicationTestUtilities::RunAdjointSensitivityDataTest<
        KEpsilonElementData::EpsilonElementData<2>,
        KEpsilonAdjointElementData::EpsilonAdjointShapeDerivatives<2, 3>::Data,
        KEpsilonAdjointElementData::EpsilonAdjointShapeDerivatives<2, 3>::ShapeDerivatives>(
        model, KEpsilonAddVariables, KEpsilonSetVariables, KEpsilonUpdateVariables, 2, 1e-6, 1e-5);
}

} // namespace Testing
} // namespace Kratos