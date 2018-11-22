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
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "processes/assign_scalar_field_to_entities_process.h"

namespace Kratos
{
template<class TEntity>
AssignScalarFieldToEntitiesProcess<TEntity>::AssignScalarFieldToEntitiesProcess(
    ModelPart& rModelPart,
    Parameters rParameters
    ) : Process() ,
        mrModelPart(rModelPart)
{
    KRATOS_TRY

    Parameters default_parameters( R"(
    {
        "model_part_name":"MODEL_PART_NAME",
        "mesh_id": 0,
        "variable_name": "VARIABLE_NAME",
        "interval"        : [0.0, 1e30],
        "value"           : "please give an expression in terms of the variable x, y, z, t",
        "local_axes" : {}
    }  )" );

    // Validate against defaults -- this ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mMeshId       = rParameters["mesh_id"].GetInt();
    mVariableName = rParameters["variable_name"].GetString();

    mpFunction = PythonGenericFunctionUtility::Pointer( new PythonGenericFunctionUtility(rParameters["value"].GetString(),  rParameters["local_axes"]));

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarFieldToEntitiesProcess<TEntity>::Execute()
{
    KRATOS_TRY;

    ProcessInfo& r_current_process_info = mrModelPart.GetProcessInfo();

    const double r_current_time = r_current_process_info[TIME];

    if( KratosComponents< Variable<double> >::Has( mVariableName ) ) { //case of scalar variable
        InternalAssignValueScalar<>(KratosComponents< Variable<double> >::Get(mVariableName), r_current_time);
    } else if( KratosComponents< array_1d_component_type >::Has( mVariableName ) ) { //case of component variable
        InternalAssignValueScalar<>(KratosComponents< array_1d_component_type >::Get(mVariableName), r_current_time);
    } else if( KratosComponents< Variable<Vector> >::Has( mVariableName ) ) { //case of vector variable
        InternalAssignValueVector<>(KratosComponents< Variable<Vector> >::Get(mVariableName), r_current_time);
    } else {
        KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mVariableName << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarFieldToEntitiesProcess<TEntity>::CallFunction(
    const typename TEntity::Pointer& pEntity,
    const double Time,
    Vector& rValue
    )
{
    GeometryType& r_entity_geometry = pEntity->GetGeometry();
    const SizeType size = r_entity_geometry.size();

    if(rValue.size() !=  size) {
        rValue.resize(size,false);
    }

    for (IndexType i=0; i<size; ++i) {
        rValue[i] = mpFunction->CallFunction(r_entity_geometry[i].X(),r_entity_geometry[i].Y(),r_entity_geometry[i].Z(),Time  );
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarFieldToEntitiesProcess<TEntity>::CallFunctionComponents(
    const typename TEntity::Pointer& pEntity,
    const double Time,
    double& rValue
    )
{
    GeometryType& r_entity_geometry = pEntity->GetGeometry();
    const array_1d<double,3>& center = r_entity_geometry.Center();

    rValue = mpFunction->CallFunction(center[0],center[1],center[2], Time  );
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarFieldToEntitiesProcess<TEntity>::CallFunctionLocalSystem(
    const typename TEntity::Pointer& pEntity,
    const double Time,
    Vector& rValue
    )
{

    GeometryType& r_entity_geometry = pEntity->GetGeometry();
    const SizeType size = r_entity_geometry.size();

    if (rValue.size() !=  size) {
        rValue.resize(size,false);
    }

    for (IndexType i=0; i<size; ++i) {
        rValue[i] = mpFunction->RotateAndCallFunction(r_entity_geometry[i].X(),r_entity_geometry[i].Y(),r_entity_geometry[i].Z(), Time  );
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarFieldToEntitiesProcess<TEntity>::CallFunctionLocalSystemComponents(
    const typename TEntity::Pointer& pEntity,
    const double Time,
    double& rValue
    )
{
    GeometryType& r_entity_geometry = pEntity->GetGeometry();

    const array_1d<double,3>& center = r_entity_geometry.Center();

    rValue = mpFunction->RotateAndCallFunction(center[0],center[1],center[2], Time  );
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarFieldToEntitiesProcess<TEntity>::AssignTimeDependentValue(
    const typename TEntity::Pointer& pEntity,
    const double Time,
    Vector& rValue,
    const double Value
    )
{
    GeometryType& r_entity_geometry = pEntity->GetGeometry();
    const SizeType size = r_entity_geometry.size();

    if(rValue.size() !=  size) {
        rValue.resize(size,false);
    }

    for(IndexType i=0; i<size; ++i) {
        rValue[i] = Value;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& AssignScalarFieldToEntitiesProcess<Condition>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh(mMeshId).Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>& AssignScalarFieldToEntitiesProcess<Element>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh(mMeshId).Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template class AssignScalarFieldToEntitiesProcess<Condition>;
template class AssignScalarFieldToEntitiesProcess<Element>;

}  // namespace Kratos.
