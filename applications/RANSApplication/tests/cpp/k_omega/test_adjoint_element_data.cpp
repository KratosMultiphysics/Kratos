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
#include "custom_elements/data_containers/k_omega/k_element_data.h"
#include "custom_elements/data_containers/k_omega/k_element_data_derivatives.h"
#include "custom_elements/data_containers/k_omega/omega_element_data.h"
#include "custom_elements/data_containers/k_omega/omega_element_data_derivatives.h"

namespace Kratos
{
namespace Testing
{

namespace
{
void KOmegaUpdateVariables(ModelPart& rModelPart)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.SetValue(TURBULENT_KINETIC_ENERGY, r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY));
        r_node.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE));
    }
}

void KOmegaAddVariables(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
    rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2);
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

void KOmegaSetVariables(ModelPart& rModelPart)
{
    using namespace RansApplicationTestUtilities;

    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY, -10.0, 10.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 1.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 1.0, 50.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, 1.0, 50.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_1, 1.0, 10.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_2, 1.0, 10.0);

    KOmegaUpdateVariables(rModelPart);

    auto& r_process_info = rModelPart.GetProcessInfo();
    r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 0.5);
    r_process_info.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA, 1.5);
    r_process_info.SetValue(TURBULENCE_RANS_BETA, 2.5);
    r_process_info.SetValue(TURBULENCE_RANS_GAMMA, 3.5);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 2.1);
}

void KOmegaSetProperties(Properties& rProperties)
{
    rProperties.SetValue(DENSITY, 1.0);
    rProperties.SetValue(DYNAMIC_VISCOSITY, 1e-2);
    rProperties.SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("RansKOmegaNewtonian2DLaw").Clone());
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(KOmegaKElementDataUDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaElementData::KElementData<2>, KOmegaElementData::KElementDataDerivatives<2, 3>::UDerivative>(
        model, KOmegaAddVariables, KOmegaSetVariables,
        KOmegaSetProperties, KOmegaUpdateVariables, 2, 1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaKElementDataKDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaElementData::KElementData<2>, KOmegaElementData::KElementDataDerivatives<2, 3>::KDerivative>(
        model, KOmegaAddVariables, KOmegaSetVariables,
        KOmegaSetProperties, KOmegaUpdateVariables, 2, 1e-6, 1e-4);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaKElementDataOmegaDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaElementData::KElementData<2>, KOmegaElementData::KElementDataDerivatives<2, 3>::OmegaDerivative>(
        model, KOmegaAddVariables, KOmegaSetVariables,
        KOmegaSetProperties, KOmegaUpdateVariables, 2, 1e-4, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaKElementDataShapeDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaElementData::KElementData<2>, KOmegaElementData::KElementDataDerivatives<2, 3>::ShapeDerivative>(
        model, KOmegaAddVariables, KOmegaSetVariables,
        KOmegaSetProperties, KOmegaUpdateVariables, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaOmegaElementDataUDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaElementData::OmegaElementData<2>, KOmegaElementData::OmegaElementDataDerivatives<2, 3>::UDerivative>(
        model, KOmegaAddVariables, KOmegaSetVariables,
        KOmegaSetProperties, KOmegaUpdateVariables, 2, 1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaOmegaElementDataKDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaElementData::OmegaElementData<2>, KOmegaElementData::OmegaElementDataDerivatives<2, 3>::KDerivative>(
        model, KOmegaAddVariables, KOmegaSetVariables,
        KOmegaSetProperties, KOmegaUpdateVariables, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaOmegaElementDataOmegaDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaElementData::OmegaElementData<2>, KOmegaElementData::OmegaElementDataDerivatives<2, 3>::OmegaDerivative>(
        model, KOmegaAddVariables, KOmegaSetVariables,
        KOmegaSetProperties, KOmegaUpdateVariables, 2, 1e-4, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaOmegaElementDataShapeDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaElementData::OmegaElementData<2>, KOmegaElementData::OmegaElementDataDerivatives<2, 3>::ShapeDerivative>(
        model, KOmegaAddVariables, KOmegaSetVariables,
        KOmegaSetProperties, KOmegaUpdateVariables, 2, 1e-6, 1e-5);
}

} // namespace Testing
} // namespace Kratos