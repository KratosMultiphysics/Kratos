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
ModelPart& CreateRansCircularConvectionRFCAdjoint2D3NModelPart(
    Model& rModel,
    const std::string& rElementName)
{
    const auto& set_variable_values = [](ModelPart& rModelPart) {
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY_POTENTIAL, 50.0, 100.0, 0);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY_POTENTIAL_RATE, 50.0, 100.0, 0);

        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY_POTENTIAL, 50.0, 100.0, 1);
        FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY_POTENTIAL_RATE, 50.0, 100.0, 1);

        auto& r_process_info = rModelPart.GetProcessInfo();
        r_process_info.SetValue(DOMAIN_SIZE, 2);
        r_process_info.SetValue(DYNAMIC_TAU, 0.04);
        r_process_info.SetValue(DELTA_TIME, 0.01);
        r_process_info.SetValue(BOSSAK_ALPHA, -0.03);
        r_process_info.SetValue(OSS_SWITCH, 0);
        r_process_info.SetValue(RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT, 1.3);
        r_process_info.SetValue(RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT, 2.8);
        r_process_info.SetValue(CIRCULAR_CONVECTION_ROTATION_CLOCKWISE, true);
        r_process_info.SetValue(CIRCULAR_CONVECTION_ROTATION_CENTER, array_1d<double, 3>({0.01, 0.01, 0.0}));
    };

    // preparing primal model part
    const auto& add_solution_step_variables = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
        rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL_RATE);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_1);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_SCALAR_1_ADJOINT_3);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUX_ADJOINT_SCALAR_1);
    };
    const auto& add_dofs = [](ModelPart::NodeType& rNode) {
        rNode.AddDof(VELOCITY_POTENTIAL);
        rNode.AddDof(RANS_SCALAR_1_ADJOINT_1);
    };

    const auto& set_element_properties = [](ModelPart::PropertiesType& rProperties) {
        Parameters cl_parameters(R"(
        {
            "name"              : "RansKEpsilonNewtonian2DLaw"
        })");
        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>::Get("RansKEpsilonNewtonian2DLaw").Create(cl_parameters, rProperties);
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
void RunRansCircularConvectionRFCAdjointTest(
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
    auto& r_primal_model_part = CreateRansCircularConvectionRFCAdjoint2D3NModelPart(model, rPrimalElementName);
    RansVariableUtilities::SetElementConstitutiveLaws(r_primal_model_part.Elements());

    // prepare adjoint model part
    auto& r_adjoint_model_part = CreateRansCircularConvectionRFCAdjoint2D3NModelPart(model, "RansCircularConvectionRFCAdjoint2D3N");

    const auto& update_function = [&](ModelPart& rModelPart) {
        const int number_of_nodes = rModelPart.NumberOfNodes();
        const double bossak_alpha = rModelPart.GetProcessInfo()[BOSSAK_ALPHA];

        for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
            auto& r_node = *(rModelPart.NodesBegin() + i_node);
            r_node.FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) = FluidAdjointTestUtilities::CalculateRelaxedVariableRate(bossak_alpha, VELOCITY_POTENTIAL_RATE, r_node);
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

KRATOS_TEST_CASE_IN_SUITE(RansCircularConvectionRFCAdjointGetDofListTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansCircularConvectionRFCAdjoint2D3NModelPart(model, "RansCircularConvectionRFCAdjoint2D3N");

    FluidTestUtilities::RunEntityGetDofListTest(
        model_part.Elements(), model_part.GetProcessInfo(),
        {&RANS_SCALAR_1_ADJOINT_1});
}

KRATOS_TEST_CASE_IN_SUITE(RansCircularConvectionRFCAdjointEquationIdVectorTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansCircularConvectionRFCAdjoint2D3NModelPart(model, "RansCircularConvectionRFCAdjoint2D3N");

    FluidTestUtilities::RunEntityEquationIdVectorTest(
        model_part.Elements(), model_part.GetProcessInfo(),
        {&RANS_SCALAR_1_ADJOINT_1});
}

KRATOS_TEST_CASE_IN_SUITE(RansCircularConvectionRFCAdjointGetValuesVectorTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansCircularConvectionRFCAdjoint2D3NModelPart(model, "RansCircularConvectionRFCAdjoint2D3N");

    FluidTestUtilities::RunEntityGetValuesVectorTest(
        model_part.Elements(),
        {&RANS_SCALAR_1_ADJOINT_1});
}

KRATOS_TEST_CASE_IN_SUITE(RansCircularConvectionRFCAdjointGetFirstDerivativesVectorTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansCircularConvectionRFCAdjoint2D3NModelPart(model, "RansCircularConvectionRFCAdjoint2D3N");

    FluidTestUtilities::RunEntityGetFirstDerivativesVectorTest(
        model_part.Elements(),
        [](const ModelPart::NodeType& rNode) -> Vector { return ZeroVector(1); });
}

KRATOS_TEST_CASE_IN_SUITE(RansCircularConvectionRFCAdjointGetSecondDerivativesVectorTest, KratosRansFastSuite)
{
    Model model;
    auto& model_part = CreateRansCircularConvectionRFCAdjoint2D3NModelPart(model, "RansCircularConvectionRFCAdjoint2D3N");

    Vector temp(1);
    FluidTestUtilities::RunEntityGetSecondDerivativesVectorTest(
        model_part.Elements(),
        [&](const ModelPart::NodeType& rNode) -> Vector {
            temp[0] = rNode.FastGetSolutionStepValue(RANS_SCALAR_1_ADJOINT_3);
            return temp;
        });
}

// /********************************************************************************************************/
// /***************************** RansCircularConvectionRFC2D3N Derivative Tests ***************************/
// /********************************************************************************************************/


KRATOS_TEST_CASE_IN_SUITE(RansCircularConvectionRFCAdjointCalculateFirstDerivativesLHS_PhiPhi, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansCircularConvectionRFCAdjointTest("RansCircularConvectionRFC2D3N", VELOCITY_POTENTIAL, derivatives_method, 0, 0, 1e-8, 1e-5);
}


KRATOS_TEST_CASE_IN_SUITE(RansCircularConvectionRFCAdjointCalculateFirstDerivativesLHS_PhiX, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansCircularConvectionRFCAdjointTest("RansCircularConvectionRFC2D3N", SHAPE_SENSITIVITY, derivatives_method, 0, 0, 1e-8, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(RansCircularConvectionRFCAdjointCalculateSecondDerivativesLHS_Phi, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ElementType& rElement,
                                        const ProcessInfo& rProcessInfo) {
        rElement.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
        const double bossak_alpha = rProcessInfo[BOSSAK_ALPHA];
        noalias(rMatrix) = rMatrix * (1.0 - bossak_alpha);
    };

    RunRansCircularConvectionRFCAdjointTest("RansCircularConvectionRFC2D3N", VELOCITY_POTENTIAL_RATE, derivatives_method, 0, 0, 1e-6, 1e-5);
}

} // namespace Testing
} // namespace Kratos