//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Josep Maria Carbonell
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "processes/assign_scalar_variable_to_entities_process.h"

namespace Kratos
{

template<class TEntity>
AssignScalarVariableToEntitiesProcess<TEntity>::AssignScalarVariableToEntitiesProcess(
    ModelPart& rModelPart,
    Parameters rParameters
    ) : Process(Flags()) ,
        mrModelPart(rModelPart)
{
    KRATOS_TRY

    Parameters default_parameters( R"(
    {
        "model_part_name":"MODEL_PART_NAME",
        "mesh_id": 0,
        "variable_name": "VARIABLE_NAME",
        "value" : 1.0
    }  )" );

    // Validate against defaults -- this ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mMeshId       = rParameters["mesh_id"].GetInt();
    mVariableName = rParameters["variable_name"].GetString();

    if( KratosComponents< Variable<double> >::Has( mVariableName )) { //case of double variable
        mDoubleValue = rParameters["value"].GetDouble();
    }
    else if( KratosComponents<array_1d_component_type>::Has( mVariableName ) ) { //case of component variable
        mDoubleValue = rParameters["value"].GetDouble();
    }
    else if( KratosComponents< Variable<int> >::Has( mVariableName ) ) { //case of int variable
        mIntValue = rParameters["value"].GetInt();
    }
    else if( KratosComponents< Variable<bool> >::Has( mVariableName ) ) { //case of bool variable
        mBoolValue = rParameters["value"].GetBool();
    } else {
        KRATOS_ERROR <<"Trying to set a variable that is not in the model_part - variable name is " << mVariableName << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarVariableToEntitiesProcess<TEntity>::Execute()
{
    KRATOS_TRY;

    if( KratosComponents< Variable<double> >::Has( mVariableName )) { //case of double variable
        InternalAssignValue<>(KratosComponents< Variable<double> >::Get(mVariableName), mDoubleValue);
    } else if( KratosComponents<array_1d_component_type>::Has( mVariableName )  ) { //case of component variable
        InternalAssignValueSerial<>(KratosComponents<array_1d_component_type>::Get(mVariableName), mDoubleValue);
    } else if( KratosComponents< Variable<int> >::Has( mVariableName ) ) { //case of int variable
        InternalAssignValue<>(KratosComponents< Variable<int> >::Get(mVariableName) , mIntValue);
    } else if( KratosComponents< Variable<bool> >::Has( mVariableName ) ) { //case of bool variable
        InternalAssignValue<>(KratosComponents< Variable<bool> >::Get(mVariableName), mBoolValue);
    } else {
        KRATOS_ERROR << "Not able to set the variable. Attempting to set variable: " << mVariableName << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& AssignScalarVariableToEntitiesProcess<Condition>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh(mMeshId).Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>& AssignScalarVariableToEntitiesProcess<Element>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh(mMeshId).Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<MasterSlaveConstraint, IndexedObject>& AssignScalarVariableToEntitiesProcess<MasterSlaveConstraint>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh(mMeshId).MasterSlaveConstraints();
}

/***********************************************************************************/
/***********************************************************************************/

template class AssignScalarVariableToEntitiesProcess<Condition>;
template class AssignScalarVariableToEntitiesProcess<Element>;
template class AssignScalarVariableToEntitiesProcess<MasterSlaveConstraint>;

}  // namespace Kratos.
