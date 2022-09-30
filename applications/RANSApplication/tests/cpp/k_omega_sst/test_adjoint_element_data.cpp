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
#include "custom_elements/data_containers/k_omega_sst/k_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/k_element_data_derivatives.h"
#include "custom_elements/data_containers/k_omega_sst/omega_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/omega_element_data_derivatives.h"

namespace Kratos
{
namespace Testing
{

namespace
{
void KOmegaSSTUpdateVariables(ModelPart& rModelPart)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.SetValue(TURBULENT_KINETIC_ENERGY, r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY));
        r_node.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE));
    }
}

void KOmegaSSTAddVariables(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(DISTANCE);
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

void KOmegaSSTSetVariables(ModelPart& rModelPart)
{
    using namespace RansApplicationTestUtilities;

    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY, -10.0, 10.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, DISTANCE, 10.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 1.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 1.0, 50.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, 1.0, 50.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_1, 1.0, 10.0);
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, RANS_AUXILIARY_VARIABLE_2, 1.0, 10.0);

    KOmegaSSTUpdateVariables(rModelPart);

    auto& r_process_info = rModelPart.GetProcessInfo();
    r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA_1, 0.5);
    r_process_info.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1, 1.5);
    r_process_info.SetValue(TURBULENCE_RANS_BETA_1, 2.5);
    r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA_2, 0.5);
    r_process_info.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2, 1.5);
    r_process_info.SetValue(TURBULENCE_RANS_BETA_2, 2.5);
    r_process_info.SetValue(TURBULENCE_RANS_A1, 3.5);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 2.1);
    r_process_info.SetValue(VON_KARMAN, 0.41);
}

void KOmegaSSTSetProperties(Properties& rProperties)
{
    rProperties.SetValue(DENSITY, 1.0);
    rProperties.SetValue(DYNAMIC_VISCOSITY, 1e-2);
    rProperties.SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("RansKOmegaSSTNewtonian2DLaw").Clone());
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTKElementDataUDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaSSTElementData::KElementData<2>, KOmegaSSTElementData::KElementDataDerivatives<2, 3>::UDerivative>(
        model, KOmegaSSTAddVariables, KOmegaSSTSetVariables,
        KOmegaSSTSetProperties, KOmegaSSTUpdateVariables, 2, 1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTKElementDataKDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaSSTElementData::KElementData<2>, KOmegaSSTElementData::KElementDataDerivatives<2, 3>::KDerivative>(
        model, KOmegaSSTAddVariables, KOmegaSSTSetVariables,
        KOmegaSSTSetProperties, KOmegaSSTUpdateVariables, 2, 1e-6, 1e-4);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTKElementDataOmegaDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaSSTElementData::KElementData<2>, KOmegaSSTElementData::KElementDataDerivatives<2, 3>::OmegaDerivative>(
        model, KOmegaSSTAddVariables, KOmegaSSTSetVariables,
        KOmegaSSTSetProperties, KOmegaSSTUpdateVariables, 2, 1e-4, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTKElementDataShapeDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaSSTElementData::KElementData<2>, KOmegaSSTElementData::KElementDataDerivatives<2, 3>::ShapeDerivative>(
        model, KOmegaSSTAddVariables, KOmegaSSTSetVariables,
        KOmegaSSTSetProperties, KOmegaSSTUpdateVariables, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTOmegaElementDataUDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaSSTElementData::OmegaElementData<2>, KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::UDerivative>(
        model, KOmegaSSTAddVariables, KOmegaSSTSetVariables,
        KOmegaSSTSetProperties, KOmegaSSTUpdateVariables, 2, 1e-7, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTOmegaElementDataKDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaSSTElementData::OmegaElementData<2>, KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::KDerivative>(
        model, KOmegaSSTAddVariables, KOmegaSSTSetVariables,
        KOmegaSSTSetProperties, KOmegaSSTUpdateVariables, 2, 1e-8, 1e-4);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTOmegaElementDataOmegaDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaSSTElementData::OmegaElementData<2>, KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::OmegaDerivative>(
        model, KOmegaSSTAddVariables, KOmegaSSTSetVariables,
        KOmegaSSTSetProperties, KOmegaSSTUpdateVariables, 2, 1e-4, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(KOmegaSSTOmegaElementDataShapeDerivative, KratosRansFastSuite)
{
    Model model;

    AdjointTestUtilities::RunElementDataDerivativeTest<
        KOmegaSSTElementData::OmegaElementData<2>, KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::ShapeDerivative>(
        model, KOmegaSSTAddVariables, KOmegaSSTSetVariables,
        KOmegaSSTSetProperties, KOmegaSSTUpdateVariables, 2, 1e-6, 1e-5);
}

} // namespace Testing
} // namespace Kratos