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

    // Validate against defaults -- this ensures no type mismatch
    const Parameters default_parameters = GetDefaultParameters();
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mMeshId       = rParameters["mesh_id"].GetInt();
    mVariableName = rParameters["variable_name"].GetString();

    mpFunction = Kratos::make_unique<GenericFunctionUtility>(rParameters["value"].GetString(), rParameters["local_axes"]);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void AssignScalarFieldToEntitiesProcess<TEntity>::Execute()
{
    KRATOS_TRY;

    const ProcessInfo& r_current_process_info = mrModelPart.GetProcessInfo();

    const double current_time = r_current_process_info[TIME];

    if( KratosComponents< Variable<double> >::Has( mVariableName ) ) { //case of scalar variable
        InternalAssignValueScalar<>(KratosComponents< Variable<double> >::Get(mVariableName), current_time);
    } else if( KratosComponents< Variable<Vector> >::Has( mVariableName ) ) { //case of vector variable
        InternalAssignValueVector<>(KratosComponents< Variable<Vector> >::Get(mVariableName), current_time);
    } else {
        KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mVariableName << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
const Parameters AssignScalarFieldToEntitiesProcess<TEntity>::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
        "model_part_name" :"MODEL_PART_NAME",
        "mesh_id"         : 0,
        "variable_name"   : "VARIABLE_NAME",
        "interval"        : [0.0, 1e30],
        "value"           : "please give an expression in terms of the variable x, y, z, t",
        "local_axes"      : {}
    }  )" );
    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node>::CallFunction(
    const typename Node::Pointer& pEntity,
    const double Time,
    Vector& rValue
    )
{
    const SizeType size = 1;

    if(rValue.size() !=  size) {
        rValue.resize(size,false);
    }

    rValue[0] = mpFunction->CallFunction(pEntity->X(),pEntity->Y(),pEntity->Z(),Time, pEntity->X0(),pEntity->Y0(),pEntity->Z0());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Condition>::CallFunction(
    const typename Condition::Pointer& pEntity,
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
        const auto& r_node = r_entity_geometry[i];
        rValue[i] = mpFunction->CallFunction(r_node.X(),r_node.Y(),r_node.Z(), Time, r_node.X0(),r_node.Y0(),r_node.Z0());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Element>::CallFunction(
    const typename Element::Pointer& pEntity,
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
        const auto& r_node = r_entity_geometry[i];
        rValue[i] = mpFunction->CallFunction(r_node.X(),r_node.Y(),r_node.Z(), Time, r_node.X0(),r_node.Y0(),r_node.Z0());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node>::CallFunctionComponents(
    const typename Node::Pointer& pEntity,
    const double Time,
    double& rValue
    )
{
    rValue = mpFunction->CallFunction(pEntity->X(),pEntity->Y(),pEntity->Z(), Time, pEntity->X0(),pEntity->Y0(),pEntity->Z0());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Condition>::CallFunctionComponents(
    const typename Condition::Pointer& pEntity,
    const double Time,
    double& rValue
    )
{
    GeometryType& r_entity_geometry = pEntity->GetGeometry();
    const array_1d<double,3>& r_center = r_entity_geometry.Center();

    rValue = mpFunction->CallFunction(r_center[0],r_center[1],r_center[2], Time);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Element>::CallFunctionComponents(
    const typename Element::Pointer& pEntity,
    const double Time,
    double& rValue
    )
{
    GeometryType& r_entity_geometry = pEntity->GetGeometry();
    const array_1d<double,3>& r_center = r_entity_geometry.Center();

    rValue = mpFunction->CallFunction(r_center[0],r_center[1],r_center[2], Time);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node>::CallFunctionLocalSystem(
    const typename Node::Pointer& pEntity,
    const double Time,
    Vector& rValue
    )
{
    const SizeType size = 1;

    if (rValue.size() !=  size) {
        rValue.resize(size,false);
    }

    rValue[0] = mpFunction->RotateAndCallFunction(pEntity->X(),pEntity->Y(),pEntity->Z(), Time, pEntity->X0(),pEntity->Y0(),pEntity->Z0());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Condition>::CallFunctionLocalSystem(
    const typename Condition::Pointer& pEntity,
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
        const auto& r_node = r_entity_geometry[i];
        rValue[i] = mpFunction->RotateAndCallFunction(r_node.X(),r_node.Y(),r_node.Z(), Time, r_node.X0(),r_node.Y0(),r_node.Z0());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Element>::CallFunctionLocalSystem(
    const typename Element::Pointer& pEntity,
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
        const auto& r_node = r_entity_geometry[i];
        rValue[i] = mpFunction->RotateAndCallFunction(r_node.X(),r_node.Y(),r_node.Z(), Time, r_node.X0(),r_node.Y0(),r_node.Z0());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node>::CallFunctionLocalSystemComponents(
    const typename Node::Pointer& pEntity,
    const double Time,
    double& rValue
    )
{
    rValue = mpFunction->RotateAndCallFunction(pEntity->X(), pEntity->Y(), pEntity->Z(), Time, pEntity->X0(), pEntity->Y0(), pEntity->Z0());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Condition>::CallFunctionLocalSystemComponents(
    const typename Condition::Pointer& pEntity,
    const double Time,
    double& rValue
    )
{
    GeometryType& r_entity_geometry = pEntity->GetGeometry();

    const array_1d<double,3>& r_center = r_entity_geometry.Center();

    rValue = mpFunction->RotateAndCallFunction(r_center[0],r_center[1],r_center[2], Time);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Element>::CallFunctionLocalSystemComponents(
    const typename Element::Pointer& pEntity,
    const double Time,
    double& rValue
    )
{
    GeometryType& r_entity_geometry = pEntity->GetGeometry();

    const array_1d<double,3>& r_center = r_entity_geometry.Center();

    rValue = mpFunction->RotateAndCallFunction(r_center[0],r_center[1],r_center[2], Time);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node>::AssignTimeDependentValue(
    const typename Node::Pointer& pEntity,
    const double Time,
    Vector& rValue,
    const double Value
    )
{
    const SizeType size = 1;

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
void AssignScalarFieldToEntitiesProcess<Condition>::AssignTimeDependentValue(
    const typename Condition::Pointer& pEntity,
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
void AssignScalarFieldToEntitiesProcess<Element>::AssignTimeDependentValue(
    const typename Element::Pointer& pEntity,
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
PointerVectorSet<Node, IndexedObject>& AssignScalarFieldToEntitiesProcess<Node>::GetEntitiesContainer()
{
    return mrModelPart.GetMesh(mMeshId).Nodes();
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

template class AssignScalarFieldToEntitiesProcess<Node>;
template class AssignScalarFieldToEntitiesProcess<Condition>;
template class AssignScalarFieldToEntitiesProcess<Element>;

}  // namespace Kratos.
