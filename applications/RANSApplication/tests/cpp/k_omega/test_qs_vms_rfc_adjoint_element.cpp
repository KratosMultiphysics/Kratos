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
#include <functional>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/constitutive_law.h"

// Application includes
#include "custom_utilities/fluid_adjoint_test_utilities.h"
#include "custom_utilities/fluid_test_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"
#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"
#include "rans_application_variables.h"

namespace Kratos
{
namespace Testing
{

namespace
{
ModelPart& CreateRansKOmegaQSVMSRFCAdjoint2D3NModelPart(
    Model& rModel,
    const std::string& rElementName)
{
    const auto& set_variable_values = [](ModelPart& rModelPart) {
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY, 50.0, 100.0, 0);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, MESH_VELOCITY, 0.0, 0.0, 0);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, PRESSURE, 5.0, 10.0, 0);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, EXTERNAL_PRESSURE, 50.0, 100.0, 0);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, ACCELERATION, 2.0, 3.0, 0);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, BODY_FORCE, 2.0, 3.0, 0);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY, 20.0, 30.0, 0);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 15.0, 25.0, 0);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 1.0, 10.0, 0);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, 50.0, 100.0, 0);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, NORMAL, 2.0, 3.0, 0);

        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY, 5.0, 10.0, 1);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, MESH_VELOCITY, 0.0, 0.0, 1);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, PRESSURE, 5.0, 10.0, 1);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, EXTERNAL_PRESSURE, 50.0, 100.0, 1);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, ACCELERATION, 2.0, 3.0, 1);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, BODY_FORCE, 2.0, 3.0, 1);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY, 2.0, 3.0, 1);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 1.0, 2.0, 1);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 10.0, 100.0, 1);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, 5.0, 10.0, 1);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, NORMAL, 2.0, 3.0, 1);

        // following values do not need to be set when OSS projections are supported by Adjoints
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, ADVPROJ, 2.0, 3.0, 0);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, DIVPROJ, 2.0, 3.0, 0);

        auto& r_process_info = rModelPart.GetProcessInfo();
        r_process_info.SetValue(DOMAIN_SIZE, 2);
        r_process_info.SetValue(DYNAMIC_TAU, 0.04);
        r_process_info.SetValue(DELTA_TIME, 0.01);
        r_process_info.SetValue(BOSSAK_ALPHA, -0.03);
        r_process_info.SetValue(OSS_SWITCH, 0);
        r_process_info.SetValue(TURBULENCE_RANS_C_MU, 3.2);
        r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 2.2);
        r_process_info.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA, 1.2);
        r_process_info.SetValue(TURBULENCE_RANS_BETA, 5.3);
        r_process_info.SetValue(TURBULENCE_RANS_GAMMA, 3.3);
        r_process_info.SetValue(RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT, 1.3);
        r_process_info.SetValue(RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT, 2.8);
    };

    // preparing primal model part
    const auto& add_solution_step_variables = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(PRESSURE);
        rModelPart.AddNodalSolutionStepVariable(ADVPROJ);
        rModelPart.AddNodalSolutionStepVariable(DIVPROJ);
        rModelPart.AddNodalSolutionStepVariable(NODAL_AREA);
        rModelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
        rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(NORMAL);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
        rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_1);
        rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_2);
        rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_3);
        rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_SCALAR_1);
        rModelPart.AddNodalSolutionStepVariable(AUX_ADJOINT_FLUID_VECTOR_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_3);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUX_ADJOINT_SCALAR_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_2_ADJOINT_3);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUX_ADJOINT_SCALAR_2);
    };
    const auto& add_dofs = [](ModelPart::NodeType& rNode) {
        rNode.AddDof(PRESSURE);
        rNode.AddDof(VELOCITY_X);
        rNode.AddDof(VELOCITY_Y);
        rNode.AddDof(VELOCITY_Z);
        rNode.AddDof(TURBULENT_KINETIC_ENERGY);
        rNode.AddDof(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);

        rNode.AddDof(ADJOINT_FLUID_SCALAR_1);
        rNode.AddDof(RANS_SCALAR_1_ADJOINT_1);
        rNode.AddDof(RANS_SCALAR_2_ADJOINT_1);
        rNode.AddDof(ADJOINT_FLUID_VECTOR_1_X);
        rNode.AddDof(ADJOINT_FLUID_VECTOR_1_Y);
        rNode.AddDof(ADJOINT_FLUID_VECTOR_1_Z);
    };

    const auto& set_element_properties = [](ModelPart::PropertiesType& rProperties) {
        Parameters cl_parameters(R"(
        {
            "name"              : "RansKOmegaNewtonian2DLaw"
        })");
        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>::Get("RansKOmegaNewtonian2DLaw").Create(cl_parameters, rProperties);
        rProperties.SetValue(CONSTITUTIVE_LAW, p_constitutive_law);
        rProperties.SetValue(DYNAMIC_VISCOSITY, 1.5);
        rProperties.SetValue(DENSITY, 1.8);
    };

    const auto& set_condition_properties = [](ModelPart::PropertiesType& rProperties) {
    };

    auto& r_model_part = FluidTestUtilities::CreateTestModelPart(
        rModel, rElementName, rElementName, "LineCondition2D2N", set_element_properties,
        set_condition_properties, add_solution_step_variables, add_dofs);
    set_variable_values(r_model_part);

    return r_model_part;
}

template<class TDataType>
void RunRansKOmegaQSVMSRFCAdjointTest(
    const std::string& rPrimalElementName,
    const Variable<TDataType>& rVariable,
    const std::function<void(Matrix&, ModelPart::ElementType&, const ProcessInfo&)>& rDerivativesRetrievalFunction,
    const IndexType EquationOffset,
    const IndexType DerivativesOffset,
    const double Delta,
    const double Tolerance)
{
    Model model;

    // prepare primal model part
    auto& r_primal_model_part = CreateRansKOmegaQSVMSRFCAdjoint2D3NModelPart(model, rPrimalElementName);
    RansVariableUtilities::SetElementConstitutiveLaws(r_primal_model_part.Elements());

    // prepare adjoint model part
    auto& r_adjoint_model_part = CreateRansKOmegaQSVMSRFCAdjoint2D3NModelPart(model, "RansKOmegaQSVMSRFCAdjoint2D3N");

    const auto& update_function = [&](ModelPart& rModelPart) {
        const int number_of_nodes = rModelPart.NumberOfNodes();
        const double bossak_alpha = rModelPart.GetProcessInfo()[BOSSAK_ALPHA];

        for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
            auto& r_node = *(rModelPart.NodesBegin() + i_node);
            r_node.SetValue(TURBULENT_KINETIC_ENERGY, r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY));
            r_node.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE));

            r_node.SetValue(RELAXED_ACCELERATION, FluidAdjointTestUtilities::CalculateRelaxedVariableRate(bossak_alpha, ACCELERATION, r_node));
            r_node.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) = FluidAdjointTestUtilities::CalculateRelaxedVariableRate(bossak_alpha, TURBULENT_KINETIC_ENERGY_RATE, r_node);
            r_node.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_2) = FluidAdjointTestUtilities::CalculateRelaxedVariableRate(bossak_alpha, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, r_node);
        }
    };

    FluidAdjointTestUtilities::RunAdjointEntityDerivativesTest(
        r_primal_model_part, r_adjoint_model_part, update_function, rVariable,
        rDerivativesRetrievalFunction, EquationOffset, DerivativesOffset, Delta, Tolerance);
}
} // namespace

/********************************************************************************************************/
/********************************************* Common Tests *********************************************/
/********************************************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointGetDofListTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansKOmegaQSVMSRFCAdjoint2D3NModelPart(model, "RansKOmegaQSVMSRFCAdjoint2D3N");

    FluidTestUtilities::RunEntityGetDofListTest(
        model_part.Elements(), model_part.GetProcessInfo(),
        {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_SCALAR_1, &RANS_SCALAR_1_ADJOINT_1, &RANS_SCALAR_2_ADJOINT_1});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointEquationIdVectorTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansKOmegaQSVMSRFCAdjoint2D3NModelPart(model, "RansKOmegaQSVMSRFCAdjoint2D3N");

    FluidTestUtilities::RunEntityEquationIdVectorTest(
        model_part.Elements(), model_part.GetProcessInfo(),
        {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_SCALAR_1, &RANS_SCALAR_1_ADJOINT_1, &RANS_SCALAR_2_ADJOINT_1});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointGetValuesVectorTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansKOmegaQSVMSRFCAdjoint2D3NModelPart(model, "RansKOmegaQSVMSRFCAdjoint2D3N");

    FluidTestUtilities::RunEntityGetValuesVectorTest(
        model_part.Elements(),
        {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_SCALAR_1, &RANS_SCALAR_1_ADJOINT_1, &RANS_SCALAR_2_ADJOINT_1});
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointGetFirstDerivativesVectorTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansKOmegaQSVMSRFCAdjoint2D3NModelPart(model, "RansKOmegaQSVMSRFCAdjoint2D3N");

    FluidTestUtilities::RunEntityGetFirstDerivativesVectorTest(
        model_part.Elements(),
        [](const ModelPart::NodeType& rNode) -> Vector { return ZeroVector(5); });
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointGetSecondDerivativesVectorTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansKOmegaQSVMSRFCAdjoint2D3NModelPart(model, "RansKOmegaQSVMSRFCAdjoint2D3N");

    Vector temp(5);
    FluidTestUtilities::RunEntityGetSecondDerivativesVectorTest(
        model_part.Elements(),
        [&](const ModelPart::NodeType& rNode) -> Vector {
            const auto& value = rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3);
            temp[0] = value[0];
            temp[1] = value[1];
            temp[2] = 0.0;
            temp[3] = rNode.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_3);
            temp[4] = rNode.FastGetSolutionStepValue(RANS_SCALAR_2_ADJOINT_3);
            return temp;
        });
}

/********************************************************************************************************/
/************************************** QSVMS2D3N Derivative Tests **************************************/
/********************************************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_UU, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("QSVMS2D3N", VELOCITY, derivatives_method, 0, 0, 1e-6, 5e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_UP, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("QSVMS2D3N", PRESSURE, derivatives_method, 0, 2, 1e-2, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_UK, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("QSVMS2D3N", TURBULENT_KINETIC_ENERGY, derivatives_method, 0, 3, 1e-5, 5e-4);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_UE, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("QSVMS2D3N", TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, derivatives_method, 0, 4, 1e-5, 1e-4);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_UX, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("QSVMS2D3N", SHAPE_SENSITIVITY, derivatives_method, 0, 0, 1e-8, 5e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateSecondDerivativesLHS_U, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rMatrix) = rMatrix * (1.0 - bossak_alpha);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("QSVMS2D3N", ACCELERATION, derivatives_method, 0, 0, 1e-2, 1e-5);
}

/********************************************************************************************************/
/********************************* RansKOmegaKRFC2D3N Derivative Tests ********************************/
/********************************************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_KU, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaKRFC2D3N", VELOCITY, derivatives_method, 3, 0, 1e-6, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_KP, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaKRFC2D3N", PRESSURE, derivatives_method, 3, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_KK, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaKRFC2D3N", TURBULENT_KINETIC_ENERGY, derivatives_method, 3, 3, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_KE, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaKRFC2D3N", TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, derivatives_method, 3, 4, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_KX, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaKRFC2D3N", SHAPE_SENSITIVITY, derivatives_method, 3, 0, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateSecondDerivativesLHS_K, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rMatrix) = rMatrix * (1.0 - bossak_alpha);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaKRFC2D3N", TURBULENT_KINETIC_ENERGY_RATE, derivatives_method, 3, 3, 1e-2, 1e-5);
}

/********************************************************************************************************/
/****************************** RansKOmegaOmegaRFC2D3N Derivative Tests *****************************/
/********************************************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_EU, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaRFC2D3N", VELOCITY, derivatives_method, 4, 0, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_EP, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaRFC2D3N", PRESSURE, derivatives_method, 4, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_EK, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaRFC2D3N", TURBULENT_KINETIC_ENERGY, derivatives_method, 4, 3, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_EE, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaRFC2D3N", TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, derivatives_method, 4, 4, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateFirstDerivativesLHS_EX, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaRFC2D3N", SHAPE_SENSITIVITY, derivatives_method, 4, 0, 1e-6, 1e-4);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaQSVMSRFCAdjointCalculateSecondDerivativesLHS_E, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rMatrix) = rMatrix * (1.0 - bossak_alpha);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaRFC2D3N", TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, derivatives_method, 4, 4, 1e-2, 3e-5);
}

} // namespace Testing
} // namespace Kratos