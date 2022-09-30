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
#include "includes/cfd_variables.h"
#include "includes/properties.h"
#include "testing/testing.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/adjoint_test_utilities.h"
#include "custom_utilities/fluid_test_utilities.h"
#include "rans_application_variables.h"

// Element data containers
#include "custom_elements/data_containers/k_epsilon/k_element_data.h"
#include "custom_elements/data_containers/k_epsilon/k_element_data_derivatives.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data_derivatives.h"

namespace Kratos
{
namespace Testing
{

namespace
{
void KEpsilonUpdateVariables(ModelPart& rModelPart)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.SetValue(TURBULENT_KINETIC_ENERGY, r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY));
        r_node.SetValue(TURBULENT_ENERGY_DISSIPATION_RATE, r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE));
    }
}

void KEpsilonAddVariables(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
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

    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY, -10.0, 10.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE, 1.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 1.0, 50.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE_2, 1.0, 50.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_1, 1.0, 10.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_2, 1.0, 10.0);

    KEpsilonUpdateVariables(rModelPart);

    auto& r_process_info = rModelPart.GetProcessInfo();
    r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 0.5);
    r_process_info.SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.5);
    r_process_info.SetValue(TURBULENCE_RANS_C1, 2.5);
    r_process_info.SetValue(TURBULENCE_RANS_C2, 3.5);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 2.1);
}

void KEpsilonSetProperties(Properties& rProperties)
{
    rProperties.SetValue(DENSITY, 1.0);
    rProperties.SetValue(DYNAMIC_VISCOSITY, 1e-2);
    rProperties.SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("RansKEpsilonNewtonian2DLaw").Clone());
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(KEpsilonKElementDataUDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KEpsilonElementData::KElementData<2>, KEpsilonElementData::KElementDataDerivatives<2, 3>::UDerivative>(
        model, KEpsilonAddVariables, KEpsilonSetVariables,
        KEpsilonSetProperties, KEpsilonUpdateVariables, 2, 1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilonKElementDataKDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KEpsilonElementData::KElementData<2>, KEpsilonElementData::KElementDataDerivatives<2, 3>::KDerivative>(
        model, KEpsilonAddVariables, KEpsilonSetVariables,
        KEpsilonSetProperties, KEpsilonUpdateVariables, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilonKElementDataEpsilonDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KEpsilonElementData::KElementData<2>, KEpsilonElementData::KElementDataDerivatives<2, 3>::EpsilonDerivative>(
        model, KEpsilonAddVariables, KEpsilonSetVariables,
        KEpsilonSetProperties, KEpsilonUpdateVariables, 2, 1e-4, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilonKElementDataShapeDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KEpsilonElementData::KElementData<2>, KEpsilonElementData::KElementDataDerivatives<2, 3>::ShapeDerivative>(
        model, KEpsilonAddVariables, KEpsilonSetVariables,
        KEpsilonSetProperties, KEpsilonUpdateVariables, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilonEpsilonElementDataUDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KEpsilonElementData::EpsilonElementData<2>, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::UDerivative>(
        model, KEpsilonAddVariables, KEpsilonSetVariables,
        KEpsilonSetProperties, KEpsilonUpdateVariables, 2, 1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilonEpsilonElementDataKDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KEpsilonElementData::EpsilonElementData<2>, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::KDerivative>(
        model, KEpsilonAddVariables, KEpsilonSetVariables,
        KEpsilonSetProperties, KEpsilonUpdateVariables, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilonEpsilonElementDataEpsilonDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KEpsilonElementData::EpsilonElementData<2>, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::EpsilonDerivative>(
        model, KEpsilonAddVariables, KEpsilonSetVariables,
        KEpsilonSetProperties, KEpsilonUpdateVariables, 2, 1e-4, 1e-4);
}

KRATOS_TEST_CASE_IN_SUITE(KEpsilonEpsilonElementDataShapeDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KEpsilonElementData::EpsilonElementData<2>, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::ShapeDerivative>(
        model, KEpsilonAddVariables, KEpsilonSetVariables,
        KEpsilonSetProperties, KEpsilonUpdateVariables, 2, 1e-6, 1e-5);
}

} // namespace Testing
} // namespace Kratos