//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "utilities/variable_utils.h"
#include "processes/apply_component_table_process.h"

namespace Kratos
{

///@name Kratos Classes
///@{
template<class TVariableType>
ApplyComponentTableProcess<TVariableType>::ApplyComponentTableProcess(
    ModelPart& rModelPart,
    Parameters rParameters)
    : Process(Flags()),
      mrModelPart(rModelPart)
{
    KRATOS_TRY

    // Only include validation with c++11 since raw_literals do not exist in c++03
    Parameters default_parameters(R"(
    {
        "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
        "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
        "time_variable_name" : "TIME",
        "is_fixed": false,
        "value" : 1.0,
        "table" : 1
    })");

    // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
    // So that an error is thrown if they don't exist
    rParameters["table"];
    rParameters["variable_name"];
    rParameters["model_part_name"];
    rParameters["time_variable_name"];

    // Now validate agains defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName = rParameters["variable_name"].GetString();
    mTimeVariableName = rParameters["time_variable_name"].GetString();
    mIsFixed = rParameters["is_fixed"].GetBool();
    mInitialValue = rParameters["value"].GetDouble();

    const unsigned int table_id = rParameters["table"].GetInt();
    mpTable = rModelPart.pGetTable(table_id);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<class TVariableType>
void ApplyComponentTableProcess<TVariableType>::Execute()
{
    ExecuteInitialize();
    ExecuteInitializeSolutionStep();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ApplyComponentTableProcess<ComponentVariableType>::ExecuteInitialize()
{
    KRATOS_TRY

    const auto& r_variable_component = KratosComponents<ComponentVariableType>::Get(mVariableName);

    ImplementationExecuteInitialize(r_variable_component);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ApplyComponentTableProcess<DoubleVariableType>::ExecuteInitialize()
{
    KRATOS_TRY

    const auto& r_variable_double = KratosComponents<DoubleVariableType>::Get(mVariableName);

    ImplementationExecuteInitialize(r_variable_double);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ApplyComponentTableProcess<ComponentVariableType>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    const auto& r_variable_component = KratosComponents<ComponentVariableType>::Get(mVariableName);

    ImplementationExecuteInitializeSolutionStep(r_variable_component);
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ApplyComponentTableProcess<DoubleVariableType>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    const auto& r_variable_double = KratosComponents<DoubleVariableType>::Get(mVariableName);

    ImplementationExecuteInitializeSolutionStep(r_variable_double);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<class TVariableType>
void ApplyComponentTableProcess<TVariableType>::ImplementationExecuteInitialize(const TVariableType& rVariable)
{
    KRATOS_TRY

    const int number_of_nodes = static_cast<int>(mrModelPart.Nodes().size());

    if (number_of_nodes != 0) {
        const auto it_node_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i<number_of_nodes; i++) {
            auto it_node = it_node_begin + i;

            if (mIsFixed) {
                it_node->Fix(rVariable);
            }

            it_node->FastGetSolutionStepValue(rVariable) = mInitialValue;
        }
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<class TVariableType>
void ApplyComponentTableProcess<TVariableType>::ImplementationExecuteInitializeSolutionStep(const TVariableType& rVariable)
{
    KRATOS_TRY

    const auto& r_time_variable = KratosComponents<DoubleVariableType>::Get(mTimeVariableName);
    const double time = mrModelPart.GetProcessInfo()[r_time_variable];
    const double value = mpTable->GetValue(time);

    const int number_of_nodes = static_cast<int>(mrModelPart.Nodes().size());

    if (number_of_nodes != 0) {
        auto& r_nodes_array = mrModelPart.Nodes();
        VariableUtils().SetScalarVar<TVariableType>(rVariable, value, r_nodes_array);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template class ApplyComponentTableProcess<ComponentVariableType>;
template class ApplyComponentTableProcess<DoubleVariableType>;

} // namespace Kratos
