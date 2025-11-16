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
#include "custom_utilities/fluid_test_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "test_utilities.h"

namespace Kratos
{
namespace KOmegaTestUtilities
{
ModelPart& RansKOmegaK2D3NSetUp(
    Model& rModel,
    const std::string& rElementName)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_1);
    };

    const auto set_properties = [](Properties& rProperties) {
        rProperties.SetValue(DENSITY, 1.0);
        rProperties.SetValue(DYNAMIC_VISCOSITY, 1e-2);
        rProperties.SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("RansKOmegaNewtonian2DLaw").Clone());
    };

    using namespace RansApplicationTestUtilities;

    auto& r_model_part = CreateScalarVariableTestModelPart(
        rModel, rElementName, "LineCondition2D2N", set_properties,
        [](Properties&) {}, add_variables_function, TURBULENT_KINETIC_ENERGY, 1);

    // set nodal historical variables
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, VELOCITY, -10.0, 10.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, MESH_VELOCITY, 10.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY_RATE, 1.0, 50.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 1.0, 1000.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, RANS_AUXILIARY_VARIABLE_1, 1.0, 10.0);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.SetValue(TURBULENT_KINETIC_ENERGY, r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY));
        r_node.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE));
    }

    // set process info variables
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(TURBULENT_KINETIC_ENERGY_SIGMA, 0.5);
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 2.1);

    return r_model_part;
}

ModelPart& RansKOmegaOmega2D3NSetUp(
    Model& rModel,
    const std::string& rElementName)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2);
        rModelPart.AddNodalSolutionStepVariable(RANS_AUXILIARY_VARIABLE_2);
    };

    const auto set_properties = [](Properties& rProperties) {
        rProperties.SetValue(DENSITY, 1.0);
        rProperties.SetValue(DYNAMIC_VISCOSITY, 1e-2);
        rProperties.SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("RansKOmegaNewtonian2DLaw").Clone());
    };

    using namespace RansApplicationTestUtilities;

    auto& r_model_part = CreateScalarVariableTestModelPart(
        rModel, rElementName, "LineCondition2D2N", set_properties, [](Properties&){}, add_variables_function,
        TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 1);

    // set nodal historical variables
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, VELOCITY, -10.0, 10.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, MESH_VELOCITY, 10.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 1.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 1.0, 1000.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, 1.0, 1000.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, RANS_AUXILIARY_VARIABLE_2, 1.0, 10.0);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.SetValue(TURBULENT_KINETIC_ENERGY, r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY));
        r_node.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node.FastGetSolutionStepValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE));
    }

    // set process info variables
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(TURBULENCE_RANS_BETA, 3.1);
    r_process_info.SetValue(TURBULENCE_RANS_GAMMA, 4.2);
    r_process_info.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA, 1.1);

    return r_model_part;
}

ModelPart& RansKOmegaOmega2D2NSetUp(
    Model& rModel,
    const std::string& rConditionName)
{
    const auto add_variables_function = [](ModelPart& rModelPart) {
        rModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(MESH_VELOCITY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_KINETIC_ENERGY);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
        rModelPart.AddNodalSolutionStepVariable(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2);
    };

    const auto set_element_properties = [](Properties& rProperties) {
        rProperties.SetValue(DENSITY, 1.0);
        rProperties.SetValue(DYNAMIC_VISCOSITY, 1e-2);
        rProperties.SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("RansKOmegaNewtonian2DLaw").Clone());
    };

    const auto set_condition_properties = [](Properties& rProperties) {
        rProperties.SetValue(WALL_SMOOTHNESS_BETA, 4.2);
        rProperties.SetValue(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT, 12.0);
    };

    using namespace RansApplicationTestUtilities;

    auto& r_model_part = CreateScalarVariableTestModelPart(
        rModel, "Element2D3N", rConditionName, set_element_properties, set_condition_properties,
        add_variables_function, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 1);

    // set nodal historical variables
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, VELOCITY, -10.0, 10.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, MESH_VELOCITY, 10.0, 100.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, TURBULENT_KINETIC_ENERGY, 10.0, 40.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, 1.0, 1000.0);
    FluidTestUtilities::RandomFillHistoricalVariable(r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, 1.0, 1000.0);

    for (auto& r_condition : r_model_part.Conditions()) {
        Vector values(2);
        values[0] = r_condition.Id();
        values[1] = 20.0 * r_condition.Id();
        r_condition.SetValue(GAUSS_RANS_Y_PLUS, values);
    }

    // set process info variables
    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info.SetValue(TURBULENCE_RANS_C_MU, 0.09);
    r_process_info.SetValue(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA, 1.1);
    r_process_info.SetValue(VON_KARMAN, 3.1);

    return r_model_part;
}
} // namespace KOmegaTestUtilities
} // namespace Kratos
