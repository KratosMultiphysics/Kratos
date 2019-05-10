//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Alejandro Cornejo
//

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

    //only include validation with c++11 since raw_literals do not exist in c++03
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
    
    unsigned int TableId = rParameters["table"].GetInt();
    mpTable = rModelPart.pGetTable(TableId);

	KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ApplyComponentTableProcess<ComponentVariableType>::ExecuteInitialize()
{
    KRATOS_TRY

	ComponentVariableType variable_component = KratosComponents<ComponentVariableType>::Get(mVariableName);
    
    const int nnodes = static_cast<int>(mrModelPart.Nodes().size());

    if (nnodes != 0) {
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i<nnodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;

            if (mIsFixed) {
                it->Fix(variable_component);
            }

            it->FastGetSolutionStepValue(variable_component) = mInitialValue;
        }
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ApplyComponentTableProcess<DoubleVariableType>::ExecuteInitialize()
{
    KRATOS_TRY
    
	DoubleVariableType variable = KratosComponents<DoubleVariableType>::Get(mVariableName);
    const int nnodes = static_cast<int>(mrModelPart.Nodes().size());

    if (nnodes != 0) {
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i < nnodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            if (mIsFixed) {
                it->Fix(variable);
            }
            it->FastGetSolutionStepValue(variable) = mInitialValue;
        }
    }
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ApplyComponentTableProcess<ComponentVariableType>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY
    
    const ComponentVariableType& variable_component = KratosComponents<ComponentVariableType>::Get(mVariableName);
	DoubleVariableType time_var_component = KratosComponents<DoubleVariableType>::Get(mTimeVariableName);
    
    const double time = mrModelPart.GetProcessInfo()[time_var_component];
    const double value = mpTable->GetValue(time);
    
    const int nnodes = static_cast<int>(mrModelPart.Nodes().size());

    if (nnodes != 0) {
        auto& r_nodes_array = mrModelPart.Nodes();
        VariableUtils().SetScalarVar<ComponentVariableType>(variable_component, value, r_nodes_array);
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
template<>
void ApplyComponentTableProcess<DoubleVariableType>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

	DoubleVariableType variable   = KratosComponents<DoubleVariableType>::Get(mVariableName);
	DoubleVariableType time_var_component = KratosComponents<DoubleVariableType>::Get(mTimeVariableName);
    
    const double time = mrModelPart.GetProcessInfo()[time_var_component];
    const double value = mpTable->GetValue(time);
    
    const int nnodes = static_cast<int>(mrModelPart.Nodes().size());

    if (nnodes != 0) {
        auto& r_nodes_array = mrModelPart.Nodes();
        VariableUtils().SetScalarVar<DoubleVariableType>(variable, value, r_nodes_array);
    }
    KRATOS_CATCH("")
}  

/***********************************************************************************/
/***********************************************************************************/

template class ApplyComponentTableProcess<ComponentVariableType>;
template class ApplyComponentTableProcess<DoubleVariableType>;

} // namespace Kratos