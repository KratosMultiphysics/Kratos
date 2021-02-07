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
#include "includes/cfd_variables.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
namespace Testing
{

namespace
{

ModelPart& CreateFluidQSVMSAdjointModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName)
{
    const auto& set_variable_values = [](ModelPart& rModelPart) {
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, VELOCITY, 50.0, 100.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, MESH_VELOCITY, 50.0, 100.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, PRESSURE, 5.0, 10.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, EXTERNAL_PRESSURE, 50.0, 100.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ACCELERATION, 2.0, 3.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, BODY_FORCE, 2.0, 3.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, VELOCITY, 5.0, 10.0, 1);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, MESH_VELOCITY, 50.0, 100.0, 1);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, PRESSURE, 5.0, 10.0, 1);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, EXTERNAL_PRESSURE, 50.0, 100.0, 1);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ACCELERATION, 2.0, 3.0, 1);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, BODY_FORCE, 2.0, 3.0, 1);

        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_SCALAR_1, 1.0, 10.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_1, 1.0, 10.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_2, 1.0, 10.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_3, 1.0, 10.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_SCALAR_1, 1.0, 10.0, 1);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_1, 1.0, 10.0, 1);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_2, 1.0, 10.0, 1);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_3, 1.0, 10.0, 1);

        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_SCALAR_1, 1.0, 10.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_1, 1.0, 10.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_2, 1.0, 10.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_3, 1.0, 10.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_SCALAR_1, 1.0, 10.0, 1);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_1, 1.0, 10.0, 1);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_2, 1.0, 10.0, 1);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADJOINT_FLUID_VECTOR_3, 1.0, 10.0, 1);

        // following values do not need to be set when OSS projections are supported by Adjoints
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, ADVPROJ, 2.0, 3.0, 0);
        FluidTestUtilities::RandomFillNodalHistoricalVariable(rModelPart, DIVPROJ, 2.0, 3.0, 0);

        auto& r_process_info = rModelPart.GetProcessInfo();
        r_process_info.SetValue(DOMAIN_SIZE, 2);
        r_process_info.SetValue(DYNAMIC_TAU, 0.04);
        r_process_info.SetValue(DELTA_TIME, 0.01);
        r_process_info.SetValue(BOSSAK_ALPHA, -0.03);
        r_process_info.SetValue(OSS_SWITCH, 0);
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
        rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_1);
        rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_2);
        rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_VECTOR_3);
        rModelPart.AddNodalSolutionStepVariable(ADJOINT_FLUID_SCALAR_1);
        rModelPart.AddNodalSolutionStepVariable(AUX_ADJOINT_FLUID_VECTOR_1);
    };
    const auto& add_dofs = [](ModelPart::NodeType& rNode) {
        rNode.AddDof(PRESSURE);
        rNode.AddDof(VELOCITY_X);
        rNode.AddDof(VELOCITY_Y);
        rNode.AddDof(VELOCITY_Z);

        rNode.AddDof(ADJOINT_FLUID_SCALAR_1);
        rNode.AddDof(ADJOINT_FLUID_VECTOR_1_X);
        rNode.AddDof(ADJOINT_FLUID_VECTOR_1_Y);
        rNode.AddDof(ADJOINT_FLUID_VECTOR_1_Z);
    };

    const auto& get_element_properties = [](ModelPart& rModelPart) {
        Properties::Pointer p_element_properties = rModelPart.CreateNewProperties(0);

        Parameters cl_parameters(R"(
        {
            "name"              : "Newtonian2DLaw"
        })");
        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>::Get("Newtonian2DLaw").Create(cl_parameters, *p_element_properties);

        p_element_properties->SetValue(CONSTITUTIVE_LAW, p_constitutive_law);
        p_element_properties->SetValue(DYNAMIC_VISCOSITY, 1.5);
        p_element_properties->SetValue(DENSITY, 1.8);

        return p_element_properties;
    };

    const auto& get_condition_properties = [](ModelPart& rModelPart) {
        Properties::Pointer p_cond_properties = rModelPart.CreateNewProperties(1);
        return p_cond_properties;
    };

    auto& r_model_part = FluidTestUtilities::CreateTestModelPart(
        rModel, rModelPartName, rElementName, "LineCondition2D2N", get_element_properties,
        get_condition_properties, add_solution_step_variables, add_dofs);
    set_variable_values(r_model_part);

    return r_model_part;
}

template<class TDataType>
void RunFluidQSVMSAdjointElementTest(
    const Variable<TDataType>& rVariable,
    const std::function<void(Matrix&, ModelPart::ElementType&, const ProcessInfo&)>& rDerivativesRetrievalFunction,
    const IndexType EquationOffset,
    const IndexType DerivativesOffset,
    const double Delta,
    const double Tolerance)
{
    Model model;

    // prepare primal model part
    auto& r_primal_model_part = CreateFluidQSVMSAdjointModelPart(model, "primal", "QSVMS2D3N");

    // prepare adjoint model part
    auto& r_adjoint_model_part = CreateFluidQSVMSAdjointModelPart(model, "adjoint", "QSVMSAdjoint2D3N");

    const auto& update_function = [](ModelPart& rModelPart) {
        const int number_of_nodes = rModelPart.NumberOfNodes();
        const double bossak_alpha = rModelPart.GetProcessInfo()[BOSSAK_ALPHA];

        for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
            auto& r_node = *(rModelPart.NodesBegin() + i_node);
            r_node.SetValue(RELAXED_ACCELERATION, FluidAdjointTestUtilities::CalculateRelaxedVariableRate(bossak_alpha, ACCELERATION, r_node));
        }
    };

    FluidAdjointTestUtilities::Testing<ModelPart::ElementsContainerType>::RunAdjointEntityDerivativesTest(
        r_primal_model_part, r_adjoint_model_part, update_function, rVariable,
        rDerivativesRetrievalFunction, EquationOffset, DerivativesOffset, Delta, Tolerance);
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(QSVMSAdjointGetDofListTest, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& model_part = CreateFluidQSVMSAdjointModelPart(model, "primal", "QSVMSAdjoint2D3N");

    FluidAdjointTestUtilities::Testing<ModelPart::ElementsContainerType>::RunAdjointEntityGetDofListTest(
        model_part,
        {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_SCALAR_1});
}

KRATOS_TEST_CASE_IN_SUITE(QSVMSAdjointEquationIdVectorTest, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& model_part = CreateFluidQSVMSAdjointModelPart(model, "primal", "QSVMSAdjoint2D3N");

    FluidAdjointTestUtilities::Testing<ModelPart::ElementsContainerType>::RunAdjointEntityEquationIdVectorTest(
        model_part,
        {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_SCALAR_1});
}

KRATOS_TEST_CASE_IN_SUITE(QSVMSAdjointGetValuesVectorTest, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& model_part = CreateFluidQSVMSAdjointModelPart(model, "primal", "QSVMSAdjoint2D3N");

    FluidAdjointTestUtilities::Testing<ModelPart::ElementsContainerType>::RunAdjointEntityGetValuesVectorTest(
        model_part,
        {&ADJOINT_FLUID_VECTOR_1_X, &ADJOINT_FLUID_VECTOR_1_Y, &ADJOINT_FLUID_SCALAR_1});
}

KRATOS_TEST_CASE_IN_SUITE(QSVMSAdjointGetFirstDerivativesVectorTest, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& model_part = CreateFluidQSVMSAdjointModelPart(model, "primal", "QSVMSAdjoint2D3N");

    FluidAdjointTestUtilities::Testing<ModelPart::ElementsContainerType>::RunAdjointEntityGetFirstDerivativesVectorTest(
        model_part,
        [](const ModelPart::NodeType& rNode) -> Vector { return ZeroVector(3); });
}

KRATOS_TEST_CASE_IN_SUITE(QSVMSAdjointGetSecondDerivativesVectorTest, FluidDynamicsApplicationFastSuite)
{
    Model model;
    auto& model_part = CreateFluidQSVMSAdjointModelPart(model, "primal", "QSVMSAdjoint2D3N");

    Vector temp(3);
    FluidAdjointTestUtilities::Testing<ModelPart::ElementsContainerType>::RunAdjointEntityGetSecondDerivativesVectorTest(
        model_part,
        [&](const ModelPart::NodeType& rNode) -> Vector {
            const auto& value = rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3);
            temp[0] = value[0];
            temp[1] = value[1];
            temp[2] = 0.0;
            return temp;
        });
}

KRATOS_TEST_CASE_IN_SUITE(QSVMSAdjointCalculateFirstDerivativesLHSVelocity, FluidDynamicsApplicationFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunFluidQSVMSAdjointElementTest(VELOCITY, derivatives_method, 0, 0, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(QSVMSAdjointCalculateFirstDerivativesLHSPressure, FluidDynamicsApplicationFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunFluidQSVMSAdjointElementTest(PRESSURE, derivatives_method, 0, 2, 1e-6, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(QSVMSAdjointCalculateSensitivityMatrixShape, FluidDynamicsApplicationFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunFluidQSVMSAdjointElementTest(SHAPE_SENSITIVITY, derivatives_method, 0, 0, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(QSVMSAdjointCalculateSecondDerivativesLHS, FluidDynamicsApplicationFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rMatrix) = rMatrix * (1.0 - bossak_alpha);
    };

    RunFluidQSVMSAdjointElementTest(ACCELERATION, derivatives_method, 0, 0, 1e-7, 1e-5);
}

} // namespace Testing
} // namespace Kratos