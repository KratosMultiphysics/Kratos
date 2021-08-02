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
#include "includes/cfd_variables.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "test_utilities.h"

namespace Kratos
{
namespace KEpsilonTestUtilities
{
ModelPart& RansKEpsilonK2D3NSetUp(
    Model& rModel,
    const std::string& rElementName)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
    };

    const auto set_properties = [](Properties& rProperties) {
        rProperties.SetValue(DENSITY, 1.0);
        rProperties.SetValue(DYNAMIC_VISCOSITY, 1e-2);
    };

    using namespace RansApplicationTestUtilities;

    auto& r_model_part = CreateScalarVariableTestModelPart(
        rModel, rElementName, "LineCondition2D2N", add_variables_function, set_properties,
        TURBULENT_KINETIC_ENERGY, 1);

    // set nodal historical variables
    RandomFillNodalHistoricalVariable(r_model_part, VELOCITY, -10.0, 10.0);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_VISCOSITY, 1e-3, 1e-1);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY_RATE, 1.0, 50.0);
    RandomFillNodalHistoricalVariable(r_model_part, RANS_AUXILIARY_VARIABLE_1, 1.0, 10.0);

    // set process info variables
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 0.5);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 2.1);

    return r_model_part;
}

ModelPart& RansKEpsilonEpsilon2D3NSetUp(
    Model& rModel,
    const std::string& rElementName)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
    };

    const auto set_properties = [](Properties& rProperties) {
        rProperties.SetValue(DENSITY, 1.0);
        rProperties.SetValue(DYNAMIC_VISCOSITY, 1e-2);
    };

    using namespace RansApplicationTestUtilities;

    auto& r_model_part = CreateScalarVariableTestModelPart(
        rModel, rElementName, "LineCondition2D2N", add_variables_function, set_properties,
        TURBULENT_ENERGY_DISSIPATION_RATE, 1);

    // set nodal historical variables
    RandomFillNodalHistoricalVariable(r_model_part, VELOCITY, -10.0, 10.0);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_VISCOSITY, 1e-3, 1e-1);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE, 1.0, 1000.0);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE_2, 1.0, 1000.0);
    RandomFillNodalHistoricalVariable(r_model_part, RANS_AUXILIARY_VARIABLE_2, 1.0, 10.0);

    // set process info variables
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(TURBULENCE_RANS_C1, 3.1);
    r_process_info.SetValue(TURBULENCE_RANS_C2, 4.2);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 1.1);
    r_process_info.SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.1);

    return r_model_part;
}

ModelPart& RansKEpsilonEpsilon2D2NSetUp(
    Model& rModel,
    const std::string& rConditionName)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_ENERGY_DISSIPATION_RATE_2);
    };

    const auto set_properties = [](Properties& rProperties) {
        rProperties.SetValue(DENSITY, 1.0);
        rProperties.SetValue(DYNAMIC_VISCOSITY, 1e-2);
        rProperties.SetValue(WALL_SMOOTHNESS_BETA, 4.2);
        rProperties.SetValue(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT, 12.0);
    };

    using namespace RansApplicationTestUtilities;

    auto& r_model_part = CreateScalarVariableTestModelPart(
        rModel, "Element2D3N", rConditionName, add_variables_function, set_properties,
        TURBULENT_ENERGY_DISSIPATION_RATE, 1);

    // set nodal historical variables
    RandomFillNodalHistoricalVariable(r_model_part, VELOCITY, -100.0, 100.0);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_VISCOSITY, 1.5, 10.2);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 10.0, 40.0);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE, 1.0, 1000.0);
    RandomFillNodalHistoricalVariable(r_model_part, TURBULENT_ENERGY_DISSIPATION_RATE_2, 1.0, 1000.0);

    RandomFillContainerVariable<ModelPart::ConditionsContainerType>(
        r_model_part, RANS_Y_PLUS, 10.0, 100.0);

    // set process info variables
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 0.09);
    r_process_info.SetValue(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, 1.1);
    r_process_info.SetValue(VON_KARMAN, 3.1);

    return r_model_part;
}
} // namespace KEpsilonTestUtilities
} // namespace Kratos
