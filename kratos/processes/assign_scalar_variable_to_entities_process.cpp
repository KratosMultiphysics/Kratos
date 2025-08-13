//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Josep Maria Carbonell
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/variable_utils.h"
#include "processes/assign_scalar_variable_to_entities_process.h"

namespace Kratos
{

template<class TEntity, bool THistorical>
AssignScalarVariableToEntitiesProcess<TEntity, THistorical>::AssignScalarVariableToEntitiesProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : AssignScalarVariableToEntitiesProcess(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()), ThisParameters)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
AssignScalarVariableToEntitiesProcess<TEntity, THistorical>::AssignScalarVariableToEntitiesProcess(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : Process(Flags()) ,
        mrModelPart(rModelPart)
{
    KRATOS_TRY

    // Validate against defaults -- this ensures no type mismatch
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName = ThisParameters["variable_name"].GetString();

    if( KratosComponents<Variable<double>>::Has( mVariableName )) { //case of double variable
        mDoubleValue = ThisParameters["value"].GetDouble();
    } else if( KratosComponents<Variable<int>>::Has( mVariableName ) ) { //case of int variable
        mIntValue = ThisParameters["value"].GetInt();
    } else if( KratosComponents<Variable<bool>>::Has( mVariableName ) ) { //case of bool variable
        mBoolValue = ThisParameters["value"].GetBool();
    } else {
        KRATOS_ERROR <<"Trying to set a variable that is not in the model_part - variable name is " << mVariableName << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
Process::Pointer AssignScalarVariableToEntitiesProcess<TEntity, THistorical>::Create(
    Model& rModel,
    Parameters ThisParameters
    )
{
    return Kratos::make_shared<AssignScalarVariableToEntitiesProcess<TEntity, THistorical>>(rModel, ThisParameters);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
void AssignScalarVariableToEntitiesProcess<TEntity, THistorical>::Execute()
{
    KRATOS_TRY;

    if(KratosComponents<Variable<double>>::Has(mVariableName)) { //case of double variable
        InternalAssignValue<>(KratosComponents<Variable<double> >::Get(mVariableName), mDoubleValue);
    } else if( KratosComponents<Variable<int>>::Has(mVariableName)) { //case of int variable
        InternalAssignValue<>(KratosComponents<Variable<int>>::Get(mVariableName) , mIntValue);
    } else if( KratosComponents<Variable<bool>>::Has(mVariableName)) { //case of bool variable
        InternalAssignValue<>(KratosComponents<Variable<bool>>::Get(mVariableName), mBoolValue);
    } else {
        KRATOS_ERROR << "Not able to set the variable. Attempting to set variable: " << mVariableName << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
const Parameters AssignScalarVariableToEntitiesProcess<TEntity, THistorical>::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
        "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "variable_name"   : "VARIABLE_NAME",
        "value"           : 1.0
    }  )" );
    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Node, IndexedObject>& AssignScalarVariableToEntitiesProcess<Node, AssignScalarVariableToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::GetEntitiesContainer()
{
    return mrModelPart.Nodes();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Node, IndexedObject>& AssignScalarVariableToEntitiesProcess<Node, AssignScalarVariableToEntitiesProcessSettings::SaveAsHistoricalVariable>::GetEntitiesContainer()
{
    return mrModelPart.Nodes();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& AssignScalarVariableToEntitiesProcess<Condition>::GetEntitiesContainer()
{
    return mrModelPart.Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>& AssignScalarVariableToEntitiesProcess<Element>::GetEntitiesContainer()
{
    return mrModelPart.Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<MasterSlaveConstraint, IndexedObject>& AssignScalarVariableToEntitiesProcess<MasterSlaveConstraint>::GetEntitiesContainer()
{
    return mrModelPart.MasterSlaveConstraints();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
template<class TVarType>
void AssignScalarVariableToEntitiesProcess<TEntity, THistorical>::InternalAssignValue(
    TVarType& rVar,
    const typename TVarType::Type Value
    )
{
    if constexpr (THistorical == AssignScalarVariableToEntitiesProcessSettings::SaveAsHistoricalVariable) {
        VariableUtils().SetVariable(rVar, Value, GetEntitiesContainer());
    } else {
        VariableUtils().SetNonHistoricalVariable(rVar, Value, GetEntitiesContainer());
    }
}

template void AssignScalarVariableToEntitiesProcess<Node, AssignScalarVariableToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::InternalAssignValue<Variable<bool>>(
    Variable<bool>& rVar,
    const bool Value
    );

template void AssignScalarVariableToEntitiesProcess<Node, AssignScalarVariableToEntitiesProcessSettings::SaveAsHistoricalVariable>::InternalAssignValue<Variable<bool>>(
    Variable<bool>& rVar,
    const bool Value
    );

template void AssignScalarVariableToEntitiesProcess<Condition>::InternalAssignValue<Variable<bool>>(
    Variable<bool>& rVar,
    const bool Value
    );

template void AssignScalarVariableToEntitiesProcess<Element>::InternalAssignValue<Variable<bool>>(
    Variable<bool>& rVar,
    const bool Value
    );

template void AssignScalarVariableToEntitiesProcess<MasterSlaveConstraint>::InternalAssignValue<Variable<bool>>(
    Variable<bool>& rVar,
    const bool Value
    );

template void AssignScalarVariableToEntitiesProcess<Node, AssignScalarVariableToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::InternalAssignValue<Variable<int>>(
    Variable<int>& rVar,
    const int Value
    );

template void AssignScalarVariableToEntitiesProcess<Node, AssignScalarVariableToEntitiesProcessSettings::SaveAsHistoricalVariable>::InternalAssignValue<Variable<int>>(
    Variable<int>& rVar,
    const int Value
    );

template void AssignScalarVariableToEntitiesProcess<Condition>::InternalAssignValue<Variable<int>>(
    Variable<int>& rVar,
    const int Value
    );

template void AssignScalarVariableToEntitiesProcess<Element>::InternalAssignValue<Variable<int>>(
    Variable<int>& rVar,
    const int Value
    );

template void AssignScalarVariableToEntitiesProcess<MasterSlaveConstraint>::InternalAssignValue<Variable<int>>(
    Variable<int>& rVar,
    const int Value
    );

template void AssignScalarVariableToEntitiesProcess<Node, AssignScalarVariableToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::InternalAssignValue<Variable<double>>(
    Variable<double>& rVar,
    const double Value
    );

template void AssignScalarVariableToEntitiesProcess<Node, AssignScalarVariableToEntitiesProcessSettings::SaveAsHistoricalVariable>::InternalAssignValue<Variable<double>>(
    Variable<double>& rVar,
    const double Value
    );

template void AssignScalarVariableToEntitiesProcess<Condition>::InternalAssignValue<Variable<double>>(
    Variable<double>& rVar,
    const double Value
    );

template void AssignScalarVariableToEntitiesProcess<Element>::InternalAssignValue<Variable<double>>(
    Variable<double>& rVar,
    const double Value
    );

template void AssignScalarVariableToEntitiesProcess<MasterSlaveConstraint>::InternalAssignValue<Variable<double>>(
    Variable<double>& rVar,
    const double Value
    );

/***********************************************************************************/
/***********************************************************************************/

template class AssignScalarVariableToEntitiesProcess<Node, AssignScalarVariableToEntitiesProcessSettings::SaveAsNonHistoricalVariable>;
template class AssignScalarVariableToEntitiesProcess<Node, AssignScalarVariableToEntitiesProcessSettings::SaveAsHistoricalVariable>;
template class AssignScalarVariableToEntitiesProcess<Condition>;
template class AssignScalarVariableToEntitiesProcess<Element>;
template class AssignScalarVariableToEntitiesProcess<MasterSlaveConstraint>;

}  // namespace Kratos.
