//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes
#include <functional>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "custom_elements/evm_k_epsilon/evm_k_epsilon_adjoint_utilities.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_logarithmic_y_plus_velocity_sensitivities_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_calculation_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_sensitivities_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_y_plus_wall_function_process.h"
#include "custom_processes/auxiliary_processes/rans_nut_y_plus_wall_function_sensitivities_process.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/test_utilities.h"
#include "test_k_epsilon_utilities.h"

namespace Kratos
{
namespace Testing
{
void InitializeModelPart(ModelPart& rModelPart)
{
    rModelPart.GetNode(1).FastGetSolutionStepValue(DENSITY) = 1.0;
    rModelPart.GetNode(2).FastGetSolutionStepValue(DENSITY) = 1.0;
    rModelPart.GetNode(3).FastGetSolutionStepValue(DENSITY) = 1.0;

    rModelPart.GetNode(1).FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 10.0;
    rModelPart.GetNode(2).FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 10.0;
    rModelPart.GetNode(3).FastGetSolutionStepValue(KINEMATIC_VISCOSITY) = 10.0;

    rModelPart.GetNode(1).FastGetSolutionStepValue(PRESSURE) = 707.3284450334143;
    rModelPart.GetNode(2).FastGetSolutionStepValue(PRESSURE) = 0.0;
    rModelPart.GetNode(3).FastGetSolutionStepValue(PRESSURE) = 1025.6378081853275;

    rModelPart.GetNode(1).FastGetSolutionStepValue(RANS_Y_PLUS) = 17.25677;
    rModelPart.GetNode(2).FastGetSolutionStepValue(RANS_Y_PLUS) = 0.0;
    rModelPart.GetNode(3).FastGetSolutionStepValue(RANS_Y_PLUS) = 0.0;

    rModelPart.GetNode(1).FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) = 0.0;
    rModelPart.GetNode(2).FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) =
        24.86820099128935;
    rModelPart.GetNode(3).FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE) =
        21.345374206136576;

    rModelPart.GetNode(1).FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 0.0;
    rModelPart.GetNode(2).FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 110.9911147569225;
    rModelPart.GetNode(3).FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY) = 0.7500000000000003;

    rModelPart.GetNode(1).FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 70.752757;
    rModelPart.GetNode(2).FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 44.58354186283689;
    rModelPart.GetNode(3).FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 0.002371708245126285;

    auto& value_1 = rModelPart.GetNode(1).FastGetSolutionStepValue(VELOCITY);
    value_1[0] = 0.0;
    value_1[1] = 0.0;
    value_1[2] = 0.0;
    auto& value_2 = rModelPart.GetNode(2).FastGetSolutionStepValue(VELOCITY);
    value_2[0] = 6.710215596607843;
    value_2[1] = 0.0;
    value_2[2] = 0.0;
    auto& value_3 = rModelPart.GetNode(3).FastGetSolutionStepValue(VELOCITY);
    value_3[0] = 10.000000000000002;
    value_3[1] = -10.000000000000002;
    value_3[2] = 0.0;

    rModelPart.GetNode(1).FastGetSolutionStepValue(VISCOSITY) = 80.752757;
    rModelPart.GetNode(2).FastGetSolutionStepValue(VISCOSITY) = 54.58354186283689;
    rModelPart.GetNode(3).FastGetSolutionStepValue(VISCOSITY) = 10.002371708245127;


    rModelPart.GetNode(1).FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE) = 0.0;
    rModelPart.GetNode(2).FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE) = 0.0;
    rModelPart.GetNode(3).FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE) = 0.0;

    rModelPart.GetNode(1).FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2) = 0.0;
    rModelPart.GetNode(2).FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2) = 0.0;
    rModelPart.GetNode(3).FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2) = 0.0;

    rModelPart.GetNode(1).FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE, 1) = 0.0;
    rModelPart.GetNode(2).FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE, 1) = 0.0;
    rModelPart.GetNode(3).FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY_RATE, 1) = 0.0;

    rModelPart.GetNode(1).FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2, 1) = 0.0;
    rModelPart.GetNode(2).FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2, 1) = 0.0;
    rModelPart.GetNode(3).FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE_2, 1) = 0.0;


    rModelPart.GetNode(1).FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) = 0.0;
    rModelPart.GetNode(2).FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) = 0.0;
    rModelPart.GetNode(3).FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_1) = 0.0;

    rModelPart.GetNode(1).FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_2) = 0.0;
    rModelPart.GetNode(2).FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_2) = 0.0;
    rModelPart.GetNode(3).FastGetSolutionStepValue(RANS_AUXILIARY_VARIABLE_2) = 0.0;

    auto& sub_model_part = rModelPart.CreateSubModelPart("boundary");
    sub_model_part.AddNode(rModelPart.pGetNode(1));
    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    r_process_info[DOMAIN_SIZE] = 2;
    r_process_info[IS_RESTARTED] = 0;
    r_process_info[BOSSAK_ALPHA] = 1.0;
    r_process_info[IS_CO_SOLVING_PROCESS_ACTIVE] = 1;
    r_process_info[DYNAMIC_TAU] = 0.0;
    r_process_info[WALL_SMOOTHNESS_BETA] = 5.2;
    r_process_info[WALL_VON_KARMAN] = 0.41;
    r_process_info[TURBULENCE_RANS_C_MU] = 0.09;
    r_process_info[TURBULENCE_RANS_C1] = 1.44;
    r_process_info[TURBULENCE_RANS_C2] = 1.92;
    r_process_info[TURBULENT_KINETIC_ENERGY_SIGMA] = 1.0;
    r_process_info[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA] = 1.3;
    r_process_info[TIME] = 1.0;
    r_process_info[DELTA_TIME] = 1.0;
    r_process_info[OSS_SWITCH] = 0;
    r_process_info[TURBULENT_VISCOSITY_MIN] = 1e-18;
    r_process_info[TURBULENT_VISCOSITY_MAX] = 100.0;
}

KRATOS_TEST_CASE_IN_SUITE(TEMPTEST_1, TEMPTEST)
{
    Model primal_model;
    ModelPart& r_primal_model_part = primal_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_primal_model_part, "RansEvmKEpsilonEpsilon2D3N");

    Model adjoint_model;
    ModelPart& r_adjoint_model_part = adjoint_model.CreateModelPart("test");
    RansEvmKEpsilonModel::GenerateRansEvmKEpsilonElementTestModelPart(
        r_adjoint_model_part, "RansEvmEpsilonAdjoint2D3N");

    InitializeModelPart(r_primal_model_part);
    InitializeModelPart(r_adjoint_model_part);

    Parameters empty_nut_parameters = Parameters(R"({
        "model_part_name" : "test"
    })");
    RansNutKEpsilonHighReSensitivitiesProcess nut_sensitivities_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess adjoint_nut_process(
        adjoint_model, empty_nut_parameters);
    RansNutKEpsilonHighReCalculationProcess primal_nut_process(primal_model, empty_nut_parameters);

    Parameters empty_nut_boundary_parameters = Parameters(R"({
        "model_part_name" : "test.boundary"
    })");
    RansNutYPlusWallFunctionSensitivitiesProcess adjoint_nut_boundary_sensitivities(
        adjoint_model, empty_nut_boundary_parameters);
    RansNutYPlusWallFunctionProcess adjoint_nut_boundary(
        adjoint_model, empty_nut_boundary_parameters);
    RansNutYPlusWallFunctionProcess primal_nut_boundary(
        primal_model, empty_nut_boundary_parameters);

    std::vector<Process*> adjoint_processes;
    adjoint_processes.push_back(&adjoint_nut_process);
    adjoint_processes.push_back(&adjoint_nut_boundary);
    adjoint_processes.push_back(&nut_sensitivities_process);
    adjoint_processes.push_back(&adjoint_nut_boundary_sensitivities);

    std::vector<Process*> primal_processes;
    primal_processes.push_back(&primal_nut_process);
    primal_processes.push_back(&primal_nut_boundary);

    auto perturb_variable = [](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE);
    };

    auto calculate_sensitivity_matrix = [](Matrix& rOutput, Element& rElement,
                                           ProcessInfo& rProcessInfo) {
        rElement.CalculateFirstDerivativesLHS(rOutput, rProcessInfo);
        KRATOS_WATCH(rOutput);
    };

    RansModellingApplicationTestUtilities::RunResidualScalarSensitivityTest(
        r_primal_model_part, r_adjoint_model_part, primal_processes,
        adjoint_processes, RansEvmKEpsilonModel::UpdateVariablesInModelPart,
        calculate_sensitivity_matrix, perturb_variable, 1e-7, 1e-5);
}
} // namespace Testing
} // namespace Kratos
