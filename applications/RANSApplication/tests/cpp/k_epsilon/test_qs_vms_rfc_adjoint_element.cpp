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
#include "custom_utilities/rans_variable_utilities.h"
#include "includes/cfd_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "rans_application_variables.h"

namespace Kratos
{
namespace Testing
{

namespace
{
ModelPart& CreateRansKEpsilonQSVMSRFCAdjoint2D3NModelPart(
    Model& rModel,
    const std::string& rElementName)
{
    const auto& set_variable_values = [](ModelPart& rModelPart) {
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, VELOCITY, 50.0, 100.0, 0);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, MESH_VELOCITY, 50.0, 100.0, 0);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, PRESSURE, 5.0, 10.0, 0);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, EXTERNAL_PRESSURE, 50.0, 100.0, 0);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ACCELERATION, 2.0, 3.0, 0);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, BODY_FORCE, 2.0, 3.0, 0);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY, 20.0, 30.0, 0);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 15.0, 25.0, 0);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE, 1.0, 10.0, 0);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE_2, 50.0, 100.0, 0);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, NORMAL, 2.0, 3.0, 0);

        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, VELOCITY, 5.0, 10.0, 1);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, MESH_VELOCITY, 50.0, 100.0, 1);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, PRESSURE, 5.0, 10.0, 1);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, EXTERNAL_PRESSURE, 50.0, 100.0, 1);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ACCELERATION, 2.0, 3.0, 1);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, BODY_FORCE, 2.0, 3.0, 1);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY, 2.0, 3.0, 1);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_KINETIC_ENERGY_RATE, 1.0, 2.0, 1);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE, 10.0, 100.0, 1);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, TURBULENT_ENERGY_DISSIPATION_RATE_2, 5.0, 10.0, 1);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, NORMAL, 2.0, 3.0, 1);

        // following values do not need to be set when OSS projections are supported by Adjoints
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADVPROJ, 2.0, 3.0, 0);
        FluidAdjointTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, DIVPROJ, 2.0, 3.0, 0);

        auto& r_process_info = rModelPart.GetProcessInfo();
        r_process_info.SetValue(DOMAIN_SIZE, 2);
        r_process_info.SetValue(DYNAMIC_TAU, 0.04);
        r_process_info.SetValue(DELTA_TIME, 0.01);
        r_process_info.SetValue(BOSSAK_ALPHA, -0.03);
        r_process_info.SetValue(OSS_SWITCH, 0);
        r_process_info.SetValue(TURBULENCE_RANS_C_MU, 3.2);
        r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 2.2);
        r_process_info.SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.2);
        r_process_info.SetValue(TURBULENCE_RANS_C1, 5.3);
        r_process_info.SetValue(TURBULENCE_RANS_C2, 3.3);
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
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2);
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
        rNode.AddDof(TURBULENT_ENERGY_DISSIPATION_RATE);

        rNode.AddDof(ADJOINT_FLUID_SCALAR_1);
        rNode.AddDof(RANS_SCALAR_1_ADJOINT_1);
        rNode.AddDof(RANS_SCALAR_2_ADJOINT_1);
        rNode.AddDof(ADJOINT_FLUID_VECTOR_1_X);
        rNode.AddDof(ADJOINT_FLUID_VECTOR_1_Y);
        rNode.AddDof(ADJOINT_FLUID_VECTOR_1_Z);
    };

    const auto& get_element_properties = [](ModelPart& rModelPart) {
        Properties::Pointer p_element_properties = rModelPart.CreateNewProperties(0);

        Parameters cl_parameters(R"(
        {
            "name"              : "RansNewtonian2DLaw"
        })");
        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>::Get("RansNewtonian2DLaw").Create(cl_parameters, *p_element_properties);

        p_element_properties->SetValue(CONSTITUTIVE_LAW, p_constitutive_law);
        p_element_properties->SetValue(DYNAMIC_VISCOSITY, 1.5);
        p_element_properties->SetValue(DENSITY, 1.8);

        return p_element_properties;
    };

    const auto& get_condition_properties = [](ModelPart& rModelPart) {
        Properties::Pointer p_cond_properties = rModelPart.CreateNewProperties(1);
        return p_cond_properties;
    };

    auto& r_model_part = FluidAdjointTestUtilities::CreateTestModelPart(
        rModel, rElementName, "LineCondition2D2N", get_element_properties,
        get_condition_properties, add_solution_step_variables, add_dofs);
    set_variable_values(r_model_part);
    return r_model_part;
}

template<class TDataType>
void RunRansKEpsilonQSVMSRFCAdjointTest(
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
    auto& r_primal_model_part = CreateRansKEpsilonQSVMSRFCAdjoint2D3NModelPart(model, rPrimalElementName);
    RansVariableUtilities::SetElementConstitutiveLaws(r_primal_model_part.Elements());

    // prepare adjoint model part
    auto& r_adjoint_model_part = CreateRansKEpsilonQSVMSRFCAdjoint2D3NModelPart(model, "RansKEpsilonQSVMSRFCAdjoint2D3N");

    const auto& update_function = [](ModelPart& rModelPart) {
        const int number_of_nodes = rModelPart.NumberOfNodes();
        const double bossak_alpha = rModelPart.GetProcessInfo()[BOSSAK_ALPHA];
        const double c_mu = rModelPart.GetProcessInfo()[TURBULENCE_RANS_C_MU];

        for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
            auto& r_node = *(rModelPart.NodesBegin() + i_node);
            r_node.SetValue(RELAXED_ACCELERATION, FluidAdjointTestUtilities::DataTypeUtilities<array_1d<double, 3>>::ComputeRelaxedVariableRate(bossak_alpha, ACCELERATION, r_node));
            r_node.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) = FluidAdjointTestUtilities::DataTypeUtilities<double>::ComputeRelaxedVariableRate(bossak_alpha, TURBULENT_KINETIC_ENERGY_RATE, r_node);
            r_node.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_2) = FluidAdjointTestUtilities::DataTypeUtilities<double>::ComputeRelaxedVariableRate(bossak_alpha, TURBULENT_ENERGY_DISSIPATION_RATE_2, r_node);

            const double k = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            const double epsilon = r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = c_mu * k * k / epsilon;
        }
    };

    FluidAdjointTestUtilities::ContainerDataTypeUtilities<ModelPart::ElementsContainerType, TDataType>::RunAdjointEntityTest(
        r_primal_model_part, r_adjoint_model_part, update_function, rVariable,
        rDerivativesRetrievalFunction, EquationOffset, DerivativesOffset, Delta, Tolerance);
}
} // namespace

/********************************************************************************************************/
/********************************************* Common Tests *********************************************/
/********************************************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointGetDofListTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansKEpsilonQSVMSRFCAdjoint2D3NModelPart(model, "RansKEpsilonQSVMSRFCAdjoint2D3N");

    FluidAdjointTestUtilities::ContainerUtilities<ModelPart::ElementsContainerType>::RunAdjointEntityGetDofListTest(
        model_part,
        {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_SCALAR_1, &RANS_SCALAR_1_ADJOINT_1, &RANS_SCALAR_2_ADJOINT_1});
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointEquationIdVectorTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansKEpsilonQSVMSRFCAdjoint2D3NModelPart(model, "RansKEpsilonQSVMSRFCAdjoint2D3N");

    FluidAdjointTestUtilities::ContainerUtilities<ModelPart::ElementsContainerType>::RunAdjointEntityEquationIdVectorTest(
        model_part,
        {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_SCALAR_1, &RANS_SCALAR_1_ADJOINT_1, &RANS_SCALAR_2_ADJOINT_1});
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointGetValuesVectorTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansKEpsilonQSVMSRFCAdjoint2D3NModelPart(model, "RansKEpsilonQSVMSRFCAdjoint2D3N");

    FluidAdjointTestUtilities::ContainerUtilities<ModelPart::ElementsContainerType>::RunAdjointEntityGetValuesVectorTest(
        model_part,
        {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_SCALAR_1, &RANS_SCALAR_1_ADJOINT_1, &RANS_SCALAR_2_ADJOINT_1});
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointGetFirstDerivativesVectorTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansKEpsilonQSVMSRFCAdjoint2D3NModelPart(model, "RansKEpsilonQSVMSRFCAdjoint2D3N");

    FluidAdjointTestUtilities::ContainerUtilities<ModelPart::ElementsContainerType>::RunAdjointEntityGetFirstDerivativesVectorTest(
        model_part,
        [](const ModelPart::NodeType& rNode) -> Vector { return ZeroVector(5); });
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointGetSecondDerivativesVectorTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansKEpsilonQSVMSRFCAdjoint2D3NModelPart(model, "RansKEpsilonQSVMSRFCAdjoint2D3N");

    Vector temp(5);
    FluidAdjointTestUtilities::ContainerUtilities<ModelPart::ElementsContainerType>::RunAdjointEntityGetSecondDerivativesVectorTest(
        model_part,
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

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_UU, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("QSVMS2D3N", VELOCITY, derivatives_method, 0, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_UP, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("QSVMS2D3N", PRESSURE, derivatives_method, 0, 2, 1e-2, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_UK, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("QSVMS2D3N", TURBULENT_KINETIC_ENERGY, derivatives_method, 0, 3, 1e-5, 3e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_UE, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("QSVMS2D3N", TURBULENT_ENERGY_DISSIPATION_RATE, derivatives_method, 0, 4, 1e-5, 3e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_UX, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("QSVMS2D3N", SHAPE_SENSITIVITY, derivatives_method, 0, 0, 1e-8, 5e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateSecondDerivativesLHS_U, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rMatrix) = rMatrix * (1.0 - bossak_alpha);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("QSVMS2D3N", ACCELERATION, derivatives_method, 0, 0, 1e-2, 1e-5);
}

/********************************************************************************************************/
/********************************* RansKEpsilonKRFC2D3N Derivative Tests ********************************/
/********************************************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_KU, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonKRFC2D3N", VELOCITY, derivatives_method, 3, 0, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_KP, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonKRFC2D3N", PRESSURE, derivatives_method, 3, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_KK, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonKRFC2D3N", TURBULENT_KINETIC_ENERGY, derivatives_method, 3, 3, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_KE, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonKRFC2D3N", TURBULENT_ENERGY_DISSIPATION_RATE, derivatives_method, 3, 4, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_KX, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonKRFC2D3N", SHAPE_SENSITIVITY, derivatives_method, 3, 0, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateSecondDerivativesLHS_K, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rMatrix) = rMatrix * (1.0 - bossak_alpha);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonKRFC2D3N", TURBULENT_KINETIC_ENERGY_RATE, derivatives_method, 3, 3, 1e-2, 1e-5);
}

/********************************************************************************************************/
/****************************** RansKEpsilonEpsilonRFC2D3N Derivative Tests *****************************/
/********************************************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_EU, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonEpsilonRFC2D3N", VELOCITY, derivatives_method, 4, 0, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_EP, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonEpsilonRFC2D3N", PRESSURE, derivatives_method, 4, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_EK, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonEpsilonRFC2D3N", TURBULENT_KINETIC_ENERGY, derivatives_method, 4, 3, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_EE, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonEpsilonRFC2D3N", TURBULENT_ENERGY_DISSIPATION_RATE, derivatives_method, 4, 4, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateFirstDerivativesLHS_EX, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonEpsilonRFC2D3N", SHAPE_SENSITIVITY, derivatives_method, 4, 0, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonQSVMSRFCAdjointCalculateSecondDerivativesLHS_E, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rMatrix) = rMatrix * (1.0 - bossak_alpha);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonEpsilonRFC2D3N", TURBULENT_ENERGY_DISSIPATION_RATE_2, derivatives_method, 4, 4, 1e-2, 3e-5);
}

} // namespace Testing
} // namespace Kratos