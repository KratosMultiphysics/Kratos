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
#include "custom_utilities/fluid_test_utilities.h"
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

void SetTkeBasedUtau(ModelPart& rModelPart) {
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY, 5.0, 10.0, 0);
}

void SetUBasedUtauLinear(ModelPart& rModelPart) {
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY, 50.0, 100.0, 0);
}

void SetUBasedUtauLogarithmic(ModelPart& rModelPart) {
    FluidTestUtilities::RandomFillHistoricalVariable(rModelPart, VELOCITY, 500.0, 1000.0, 0);
}

ModelPart& CreateRansKOmegaVMSKBasedOmegaKBased2D2NModelPart(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rConditionName,
    const std::function<void(ModelPart&)>& rVelocitySetMethod)
{
    const auto& set_variable_values = [&](ModelPart& rModelPart) {
        rVelocitySetMethod(rModelPart);
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
        r_process_info.SetValue(DELTA_TIME, 0.01);
        r_process_info.SetValue(TURBULENCE_RANS_C_MU, 3.2);
        r_process_info.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA, 1.2);
        r_process_info.SetValue(VON_KARMAN, 0.41);

        VariableUtils().SetHistoricalVariableToZero(EXTERNAL_PRESSURE, rModelPart.Nodes());
        VariableUtils().SetNonHistoricalVariable(GAUSS_RANS_Y_PLUS, Vector(2), rModelPart.Conditions());
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

    const auto& set_element_properties = [](Properties& rProperties) {
        Parameters cl_parameters(R"(
        {
            "name"              : "RansKOmegaNewtonian2DLaw"
        })");
        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>::Get("RansKOmegaNewtonian2DLaw").Create(cl_parameters, rProperties);

        rProperties.SetValue(CONSTITUTIVE_LAW, p_constitutive_law);
        rProperties.SetValue(DYNAMIC_VISCOSITY, 1.5);
        rProperties.SetValue(DENSITY, 1.8);
    };

    const auto& set_condition_properties = [](Properties& rProperties) {
        rProperties.SetValue(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT, RansCalculationUtilities::CalculateLogarithmicYPlusLimit(0.41, 5.2));
        rProperties.SetValue(WALL_SMOOTHNESS_BETA, 5.2);
    };

    auto& r_model_part = FluidTestUtilities::CreateTestModelPart(
        rModel, rModelPartName, "RansKOmegaQSVMSRFCAdjoint2D3N", rConditionName, set_element_properties,
        set_condition_properties, add_solution_step_variables, add_dofs);
    set_variable_values(r_model_part);

    r_model_part.GetCondition(1).SetValue(RANS_IS_WALL_FUNCTION_ACTIVE, 1);
    r_model_part.GetCondition(1).Set(SLIP, true);

    NormalCalculationUtils().CalculateOnSimplex(r_model_part.Conditions(), 2);

    return r_model_part;
}

template<class TDataType>
void RunRansKOmegaQSVMSRFCAdjointTest(
    const std::string& rPrimalConditionName,
    const Variable<TDataType>& rVariable,
    const std::function<void(Matrix&, ModelPart::ConditionType&, const ProcessInfo&)>& rDerivativesRetrievalFunction,
    const std::function<void(ModelPart&)>& rVelocitySettingMethod,
    const IndexType EquationOffset,
    const IndexType DerivativesOffset,
    const double Delta,
    const double Tolerance)
{
    Model model;

    // prepare primal model part
    auto& r_primal_model_part = CreateRansKOmegaVMSKBasedOmegaKBased2D2NModelPart(model, "primal", rPrimalConditionName, rVelocitySettingMethod);
    RansVariableUtilities::SetElementConstitutiveLaws(r_primal_model_part.Elements());

    // prepare adjoint model part
    auto& r_adjoint_model_part = CreateRansKOmegaVMSKBasedOmegaKBased2D2NModelPart(model, "adjoint", "RansKOmegaVMSKBasedOmegaKBasedWallAdjoint2D2N", rVelocitySettingMethod);
    r_adjoint_model_part.GetElement(1).Initialize(r_adjoint_model_part.GetProcessInfo());
    NormalCalculationUtils().CalculateNormalShapeDerivativesOnSimplex(r_adjoint_model_part.Conditions(), 2);

    const auto& update_function = [&](ModelPart& rModelPart) {
        auto& r_process_info = rModelPart.GetProcessInfo();
        const double kappa = r_process_info[VON_KARMAN];

        NormalCalculationUtils().CalculateOnSimplex(rModelPart.Conditions(), 2);

        Vector Ws;
        Matrix Ns;

        const int number_of_conditions = rModelPart.NumberOfConditions();
        for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond) {
            auto& r_cond = *(rModelPart.ConditionsBegin() + i_cond);
            if (RansCalculationUtilities::IsWallFunctionActive(r_cond.GetGeometry())) {
                const auto& r_parent_element = r_cond.GetValue(NEIGHBOUR_ELEMENTS)[0];
                const auto& r_element_properties = r_parent_element.GetProperties();
                const double rho = r_element_properties[DENSITY];
                const double nu = r_element_properties[DYNAMIC_VISCOSITY] / rho;

                const auto& r_condition_properties = r_cond.GetProperties();
                const double beta = r_condition_properties[WALL_SMOOTHNESS_BETA];
                const double y_plus_limit = r_condition_properties[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];

                array_1d<double, 3> normal = r_cond.GetValue(NORMAL);
                normal /= norm_2(normal);

                const auto& r_condition_geometry = r_cond.GetGeometry();
                const double y = inner_prod(r_condition_geometry.Center() - r_parent_element.GetGeometry().Center(), normal);
                r_cond.SetValue(DISTANCE, y);

                RansCalculationUtilities::CalculateConditionGeometryData(
                    r_condition_geometry, GeometryData::IntegrationMethod::GI_GAUSS_2, Ws, Ns);
                const int number_of_gauss_points = Ws.size();

                Vector& r_gauss_y_plus = r_cond.GetValue(GAUSS_RANS_Y_PLUS);

                array_1d<double, 3> wall_velocity;
                for (int g = 0; g < number_of_gauss_points; ++g) {
                    double& y_plus = r_gauss_y_plus[g];
                    const Vector& N = row(Ns, g);

                    FluidCalculationUtilities::EvaluateInPoint(r_condition_geometry, N,
                        std::tie(wall_velocity, VELOCITY));

                    const double wall_velocity_magnitude = norm_2(wall_velocity);
                    double u_tau;
                    RansCalculationUtilities::CalculateYPlusAndUtau(y_plus, u_tau, wall_velocity_magnitude, y, nu, kappa, beta);
                    y_plus = std::max(y_plus, y_plus_limit);
                }
            }
        }

        if (rVariable == SHAPE_SENSITIVITY) {
            r_process_info.SetValue(FRACTIONAL_STEP, 200);
        }
    };

    FluidAdjointTestUtilities::RunAdjointEntityDerivativesTest(
        r_primal_model_part, r_adjoint_model_part, update_function, rVariable,
        rDerivativesRetrievalFunction, EquationOffset, DerivativesOffset, Delta, Tolerance);
}
} // namespace

/********************************************************************************************************/
/************** RansKOmegaVMSMonolithicKBasedOmegaKBasedWallCondition Derivative Tests **************/
/********************************************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UU_1, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", VELOCITY, derivatives_method, SetUBasedUtauLinear, 0, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UU_2, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", VELOCITY, derivatives_method, SetUBasedUtauLogarithmic, 0, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UU_3, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", VELOCITY, derivatives_method, SetTkeBasedUtau, 0, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UP_1, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", PRESSURE, derivatives_method, SetUBasedUtauLinear, 0, 2, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UP_2, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", PRESSURE, derivatives_method, SetUBasedUtauLogarithmic, 0, 2, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UP_3, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", PRESSURE, derivatives_method, SetTkeBasedUtau, 0, 2, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UK_1, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", TURBULENT_KINETIC_ENERGY, derivatives_method, SetUBasedUtauLinear, 0, 3, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UK_2, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", TURBULENT_KINETIC_ENERGY, derivatives_method, SetUBasedUtauLogarithmic, 0, 3, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UK_3, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", TURBULENT_KINETIC_ENERGY, derivatives_method, SetTkeBasedUtau, 0, 3, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UE_1, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, derivatives_method, SetUBasedUtauLinear, 0, 4, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UE_2, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, derivatives_method, SetUBasedUtauLogarithmic, 0, 4, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UE_3, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, derivatives_method, SetTkeBasedUtau, 0, 4, 1e-5, 1e-5);
}


KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UX_1, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", SHAPE_SENSITIVITY, derivatives_method, SetUBasedUtauLinear, 0, 0, 1e-10, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UX_2, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", SHAPE_SENSITIVITY, derivatives_method, SetUBasedUtauLogarithmic, 0, 0, 1e-10, 1e-4);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_UX_3, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", SHAPE_SENSITIVITY, derivatives_method, SetTkeBasedUtau, 0, 0, 1e-10, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionSecondDerivativesLHS_U_1, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", ACCELERATION, derivatives_method, SetUBasedUtauLinear, 0, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionSecondDerivativesLHS_U_2, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", ACCELERATION, derivatives_method, SetUBasedUtauLogarithmic, 0, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionSecondDerivativesLHS_U_3, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateSecondDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansVMSMonolithicKBasedWall2D2N", ACCELERATION, derivatives_method, SetTkeBasedUtau, 0, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EU_1, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", VELOCITY, derivatives_method, SetUBasedUtauLinear, 4, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EU_2, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", VELOCITY, derivatives_method, SetUBasedUtauLogarithmic, 4, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EU_3, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", VELOCITY, derivatives_method, SetTkeBasedUtau, 4, 0, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EP_1, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", PRESSURE, derivatives_method, SetUBasedUtauLinear, 4, 2, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EP_2, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", PRESSURE, derivatives_method, SetUBasedUtauLogarithmic, 4, 2, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EP_3, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", PRESSURE, derivatives_method, SetTkeBasedUtau, 4, 2, 1e-5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EK_1, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", TURBULENT_KINETIC_ENERGY, derivatives_method, SetUBasedUtauLinear, 4, 3, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EK_2, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", TURBULENT_KINETIC_ENERGY, derivatives_method, SetUBasedUtauLogarithmic, 4, 3, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EK_3, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", TURBULENT_KINETIC_ENERGY, derivatives_method, SetTkeBasedUtau, 4, 3, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EE_1, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, derivatives_method, SetUBasedUtauLinear, 4, 4, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EE_2, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, derivatives_method, SetUBasedUtauLogarithmic, 4, 4, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EE_3, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateFirstDerivativesLHS(rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, derivatives_method, SetTkeBasedUtau, 4, 4, 1e-7, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EX_1, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", SHAPE_SENSITIVITY, derivatives_method, SetUBasedUtauLinear, 4, 0, 1e-10, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EX_2, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", SHAPE_SENSITIVITY, derivatives_method, SetUBasedUtauLogarithmic, 4, 0, 1e-10, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaVMSMonolithicKBasedOmegaKBasedWallConditionFirstDerivativesLHS_EX_3, KratosRansFastSuite)
{
    const auto& derivatives_method = [](Matrix& rMatrix, ModelPart::ConditionType& rCondition,
                                        const ProcessInfo& rProcessInfo) {
        rCondition.CalculateSensitivityMatrix(SHAPE_SENSITIVITY, rMatrix, rProcessInfo);
    };

    RunRansKOmegaQSVMSRFCAdjointTest("RansKOmegaOmegaKBasedWall2D2N", SHAPE_SENSITIVITY, derivatives_method, SetTkeBasedUtau, 4, 0, 1e-10, 1e-5);
}

} // namespace Testing
} // namespace Kratos