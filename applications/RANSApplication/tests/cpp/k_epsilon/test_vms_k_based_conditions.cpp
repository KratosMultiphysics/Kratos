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
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/fluid_adjoint_test_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "includes/cfd_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "rans_application_variables.h"
#include "utilities/normal_calculation_utils.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
namespace Testing
{

namespace
{
ModelPart& CreateRansKEpsilonVMSKBasedEpsilonKBased2D2NModelPart(
    Model& rModel,
    const std::string& rConditionName)
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
        r_process_info.SetValue(DELTA_TIME, 0.01);
        r_process_info.SetValue(TURBULENCE_RANS_C_MU, 3.2);
        r_process_info.SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.2);
        r_process_info.SetValue(VON_KARMAN, 0.41);

        VariableUtils().SetHistoricalVariableToZero(EXTERNAL_PRESSURE, rModelPart.Nodes());
    };

    // preparing primal model part
    const auto& add_solution_step_variables = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(PRESSURE);
        rModelPart.AddNodalSolutionStepVariable(ADVPROJ);
        rModelPart.AddNodalSolutionStepVariable(DIVPROJ);
        rModelPart.AddNodalSolutionStepVariable(DENSITY);
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
        p_cond_properties->SetValue(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT, RansCalculationUtilities::CalculateLogarithmicYPlusLimit(0.41, 5.2));
        p_cond_properties->SetValue(WALL_SMOOTHNESS_BETA, 5.2);
        return p_cond_properties;
    };

    auto& r_model_part = FluidAdjointTestUtilities::CreateTestModelPart(
        rModel, "RansKEpsilonQSVMSRFCAdjoint2D3N", rConditionName, get_element_properties,
        get_condition_properties, add_solution_step_variables, add_dofs);
    set_variable_values(r_model_part);

    r_model_part.GetCondition(1).SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_model_part.GetCondition(1).Set(SLIP, true);

    NormalCalculationUtils().CalculateOnSimplex(r_model_part.Conditions(), 2);

    return r_model_part;
}

template<class TDataType>
void RunRansKEpsilonQSVMSRFCAdjointTest(
    const std::string& rPrimalConditionName,
    const Variable<TDataType>& rVariable,
    const std::function<void(Matrix&, ModelPart::ConditionType&, const ProcessInfo&)>& rDerivativesRetrievalFunction,
    const IndexType EquationOffset,
    const IndexType DerivativesOffset,
    const double Delta,
    const double Tolerance)
{
    Model model;

    // prepare primal model part
    auto& r_primal_model_part = CreateRansKEpsilonVMSKBasedEpsilonKBased2D2NModelPart(model, rPrimalConditionName);
    RansVariableUtilities::SetElementConstitutiveLaws(r_primal_model_part.Elements());

    // prepare adjoint model part
    auto& r_adjoint_model_part = CreateRansKEpsilonVMSKBasedEpsilonKBased2D2NModelPart(model, "RansKEpsilonVMSKBasedEpsilonKBasedWallAdjoint2D2N");
    r_adjoint_model_part.GetElement(1).Initialize(r_adjoint_model_part.GetProcessInfo());
    NormalCalculationUtils().CalculateNormalShapeDerivativesOnSimplex(r_adjoint_model_part.Conditions(), 2);

    const auto& update_function = [&](ModelPart& rModelPart) {
        const int number_of_nodes = rModelPart.NumberOfNodes();

        auto& r_process_info = rModelPart.GetProcessInfo();
        const double c_mu = r_process_info[TURBULENCE_RANS_C_MU];
        const double kappa = r_process_info[VON_KARMAN];

        for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
            auto& r_node = *(rModelPart.NodesBegin() + i_node);
            const double k = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY);
            const double epsilon = r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
            r_node.FastGetSolutionStepValue(TURBULENT_VISCOSITY) = c_mu * k * k / epsilon;
        }

        NormalCalculationUtils().CalculateOnSimplex(rModelPart.Conditions(), 2);

        const int number_of_conditions = rModelPart.NumberOfConditions();
        Vector N(2);
        N[0] = 0.5;
        N[1] = 0.5;
        const double nu = 1.5 / 1.8;
        for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond) {
            auto& r_cond = *(rModelPart.ConditionsBegin() + i_cond);
            if (RansCalculationUtilities::IsWallFunctionActive(r_cond)) {
                const double beta = r_cond.GetProperties()[WALL_SMOOTHNESS_BETA];

                array_1d<double, 3> wall_velocity;
                FluidCalculationUtilities::EvaluateInPoint(
                    r_cond.GetGeometry(), N, std::tie(wall_velocity, VELOCITY));

                const double wall_velocity_magnitude = norm_2(wall_velocity);
                array_1d<double, 3> normal = r_cond.GetValue(NORMAL);
                normal /= norm_2(normal);

                const double y = inner_prod(
                    r_cond.GetGeometry().Center() -
                        r_cond.GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry().Center(),
                    normal);

                double y_plus, u_tau;
                RansCalculationUtilities::CalculateYPlusAndUtau(
                    y_plus, u_tau, wall_velocity_magnitude, y, nu, kappa, beta);

                r_cond.SetValue(RANS_Y_PLUS, y_plus);
                r_cond.SetValue(FRICTION_VELOCITY, wall_velocity * u_tau / wall_velocity_magnitude);
            }
        }

        if (rVariable == SHAPE_SENSITIVITY) {
            r_process_info.SetValue(FRACTIONAL_STEP, 200);
        }
    };

    FluidAdjointTestUtilities::ContainerDataTypeUtilities<ModelPart::ConditionsContainerType, TDataType>::RunAdjointEntityTest(
        r_primal_model_part, r_adjoint_model_part, update_function, rVariable,
        rDerivativesRetrievalFunction, EquationOffset, DerivativesOffset, Delta, Tolerance);
}
} // namespace

/********************************************************************************************************/
/************** RansKEpsilonVMSMonolithicKBasedEpsilonKBasedWallCondition Derivative Tests **************/
/********************************************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonVMSMonolithicKBasedEpsilonKBasedWallConditionFirstDerivativesLHS_UU, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", VELOCITY, derivatives_method, 0, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonVMSMonolithicKBasedEpsilonKBasedWallConditionFirstDerivativesLHS_UP, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", PRESSURE, derivatives_method, 0, 2, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonVMSMonolithicKBasedEpsilonKBasedWallConditionFirstDerivativesLHS_UK, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", TURBULENT_KINETIC_ENERGY, derivatives_method, 0, 3, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonVMSMonolithicKBasedEpsilonKBasedWallConditionFirstDerivativesLHS_UE, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", TURBULENT_ENERGY_DISSIPATION_RATE, derivatives_method, 0, 4, 1e-5, 1e-5);
}


KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonVMSMonolithicKBasedEpsilonKBasedWallConditionFirstDerivativesLHS_UX, KratosRansFastSuite1)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
        KRATOS_WATCH(rMatrix);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", SHAPE_SENSITIVITY, derivatives_method, 0, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonVMSMonolithicKBasedEpsilonKBasedWallConditionSecondDerivativesLHS_U, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", ACCELERATION, derivatives_method, 0, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonVMSMonolithicKBasedEpsilonKBasedWallConditionFirstDerivativesLHS_EU, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonEpsilonKBasedWall2D2N", VELOCITY, derivatives_method, 4, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonVMSMonolithicKBasedEpsilonKBasedWallConditionFirstDerivativesLHS_EP, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonEpsilonKBasedWall2D2N", PRESSURE, derivatives_method, 4, 2, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonVMSMonolithicKBasedEpsilonKBasedWallConditionFirstDerivativesLHS_EK, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonEpsilonKBasedWall2D2N", TURBULENT_KINETIC_ENERGY, derivatives_method, 4, 3, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKEpsilonVMSMonolithicKBasedEpsilonKBasedWallConditionFirstDerivativesLHS_EE, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKEpsilonQSVMSRFCAdjointTest("RansKEpsilonEpsilonKBasedWall2D2N", TURBULENT_ENERGY_DISSIPATION_RATE, derivatives_method, 4, 4, 1e-7, 1e-5);
}

} // namespace Testing
} // namespace Kratos