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
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/assign_scalar_field_to_entities_process.h"

namespace Kratos
{
template<class TEntity, bool THistorical>
AssignScalarFieldToEntitiesProcess<TEntity, THistorical>::AssignScalarFieldToEntitiesProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : AssignScalarFieldToEntitiesProcess(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()), ThisParameters)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
AssignScalarFieldToEntitiesProcess<TEntity, THistorical>::AssignScalarFieldToEntitiesProcess(
    ModelPart& rModelPart,
    Parameters rParameters
    ) : Process() ,
        mrModelPart(rModelPart)
{
    KRATOS_TRY

    // Validate against defaults -- this ensures no type mismatch
    const Parameters default_parameters = GetDefaultParameters();
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName = rParameters["variable_name"].GetString();

    mpFunction = Kratos::make_unique<GenericFunctionUtility>(rParameters["value"].GetString(), rParameters["local_axes"]);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
Process::Pointer AssignScalarFieldToEntitiesProcess<TEntity, THistorical>::Create(
    Model& rModel,
    Parameters ThisParameters)
{
    return Kratos::make_shared<AssignScalarFieldToEntitiesProcess<TEntity, THistorical>>(rModel, ThisParameters);

}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
void AssignScalarFieldToEntitiesProcess<TEntity, THistorical>::Execute()
{
    KRATOS_TRY;

    const ProcessInfo& r_current_process_info = mrModelPart.GetProcessInfo();

    const double current_time = r_current_process_info[TIME];

    if( KratosComponents<Variable<double>>::Has(mVariableName)) { //case of scalar variable
        InternalAssignValueScalar(KratosComponents<Variable<double>>::Get(mVariableName), current_time);
    } else if( KratosComponents<Variable<Vector>>::Has(mVariableName)) { //case of vector variable
        InternalAssignValueVector(KratosComponents<Variable<Vector>>::Get(mVariableName), current_time);
    } else {
        KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mVariableName << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
const Parameters AssignScalarFieldToEntitiesProcess<TEntity, THistorical>::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
        "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
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
void AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::CallFunction(
    const typename Node::Pointer pEntity,
    const double Time,
    Vector& rValue,
    GenericFunctionUtility& rFunction
    )
{
    if(rValue.size() != 1) {
        rValue.resize(1,false);
    }

    rValue[0] = rFunction.CallFunction(pEntity->X(),pEntity->Y(),pEntity->Z(),Time, pEntity->X0(),pEntity->Y0(),pEntity->Z0());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsHistoricalVariable>::CallFunction(
    const typename Node::Pointer pEntity,
    const double Time,
    Vector& rValue,
    GenericFunctionUtility& rFunction
    )
{
    if(rValue.size() != 1) {
        rValue.resize(1,false);
    }

    rValue[0] = rFunction.CallFunction(pEntity->X(),pEntity->Y(),pEntity->Z(),Time, pEntity->X0(),pEntity->Y0(),pEntity->Z0());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Condition>::CallFunction(
    const typename Condition::Pointer pEntity,
    const double Time,
    Vector& rValue,
    GenericFunctionUtility& rFunction
    )
{
    auto& r_entity_geometry = pEntity->GetGeometry();
    const SizeType size = r_entity_geometry.size();

    if(rValue.size() !=  size) {
        rValue.resize(size,false);
    }

    for (IndexType i=0; i<size; ++i) {
        const auto& r_node = r_entity_geometry[i];
        rValue[i] = rFunction.CallFunction(r_node.X(),r_node.Y(),r_node.Z(), Time, r_node.X0(),r_node.Y0(),r_node.Z0());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Element>::CallFunction(
    const typename Element::Pointer pEntity,
    const double Time,
    Vector& rValue,
    GenericFunctionUtility& rFunction
    )
{
    auto& r_entity_geometry = pEntity->GetGeometry();
    const SizeType size = r_entity_geometry.size();

    if(rValue.size() !=  size) {
        rValue.resize(size,false);
    }

    for (IndexType i=0; i<size; ++i) {
        const auto& r_node = r_entity_geometry[i];
        rValue[i] = rFunction.CallFunction(r_node.X(),r_node.Y(),r_node.Z(), Time, r_node.X0(),r_node.Y0(),r_node.Z0());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::CallFunctionComponents(
    const typename Node::Pointer pEntity,
    const double Time,
    double& rValue,
    GenericFunctionUtility& rFunction
    )
{
    rValue = rFunction.CallFunction(pEntity->X(),pEntity->Y(),pEntity->Z(), Time, pEntity->X0(),pEntity->Y0(),pEntity->Z0());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsHistoricalVariable>::CallFunctionComponents(
    const typename Node::Pointer pEntity,
    const double Time,
    double& rValue,
    GenericFunctionUtility& rFunction
    )
{
    rValue = rFunction.CallFunction(pEntity->X(),pEntity->Y(),pEntity->Z(), Time, pEntity->X0(),pEntity->Y0(),pEntity->Z0());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Condition>::CallFunctionComponents(
    const typename Condition::Pointer pEntity,
    const double Time,
    double& rValue,
    GenericFunctionUtility& rFunction
    )
{
    auto& r_entity_geometry = pEntity->GetGeometry();
    const array_1d<double,3>& r_center = r_entity_geometry.Center();

    rValue = rFunction.CallFunction(r_center[0],r_center[1],r_center[2], Time);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Element>::CallFunctionComponents(
    const typename Element::Pointer pEntity,
    const double Time,
    double& rValue,
    GenericFunctionUtility& rFunction
    )
{
    auto& r_entity_geometry = pEntity->GetGeometry();
    const array_1d<double,3>& r_center = r_entity_geometry.Center();

    rValue = rFunction.CallFunction(r_center[0],r_center[1],r_center[2], Time);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::CallFunctionLocalSystem(
    const typename Node::Pointer pEntity,
    const double Time,
    Vector& rValue,
    GenericFunctionUtility& rFunction
    )
{
    if (rValue.size() !=  1) {
        rValue.resize(1,false);
    }

    rValue[0] = rFunction.RotateAndCallFunction(pEntity->X(),pEntity->Y(),pEntity->Z(), Time, pEntity->X0(),pEntity->Y0(),pEntity->Z0());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsHistoricalVariable>::CallFunctionLocalSystem(
    const typename Node::Pointer pEntity,
    const double Time,
    Vector& rValue,
    GenericFunctionUtility& rFunction
    )
{
    if (rValue.size() !=  1) {
        rValue.resize(1,false);
    }

    rValue[0] = rFunction.RotateAndCallFunction(pEntity->X(),pEntity->Y(),pEntity->Z(), Time, pEntity->X0(),pEntity->Y0(),pEntity->Z0());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Condition>::CallFunctionLocalSystem(
    const typename Condition::Pointer pEntity,
    const double Time,
    Vector& rValue,
    GenericFunctionUtility& rFunction
    )
{

    auto& r_entity_geometry = pEntity->GetGeometry();
    const SizeType size = r_entity_geometry.size();

    if (rValue.size() !=  size) {
        rValue.resize(size,false);
    }

    for (IndexType i=0; i<size; ++i) {
        const auto& r_node = r_entity_geometry[i];
        rValue[i] = rFunction.RotateAndCallFunction(r_node.X(),r_node.Y(),r_node.Z(), Time, r_node.X0(),r_node.Y0(),r_node.Z0());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Element>::CallFunctionLocalSystem(
    const typename Element::Pointer pEntity,
    const double Time,
    Vector& rValue,
    GenericFunctionUtility& rFunction
    )
{

    auto& r_entity_geometry = pEntity->GetGeometry();
    const SizeType size = r_entity_geometry.size();

    if (rValue.size() !=  size) {
        rValue.resize(size,false);
    }

    for (IndexType i=0; i<size; ++i) {
        const auto& r_node = r_entity_geometry[i];
        rValue[i] = rFunction.RotateAndCallFunction(r_node.X(),r_node.Y(),r_node.Z(), Time, r_node.X0(),r_node.Y0(),r_node.Z0());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::CallFunctionLocalSystemComponents(
    const typename Node::Pointer pEntity,
    const double Time,
    double& rValue,
    GenericFunctionUtility& rFunction
    )
{
    rValue = rFunction.RotateAndCallFunction(pEntity->X(), pEntity->Y(), pEntity->Z(), Time, pEntity->X0(), pEntity->Y0(), pEntity->Z0());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsHistoricalVariable>::CallFunctionLocalSystemComponents(
    const typename Node::Pointer pEntity,
    const double Time,
    double& rValue,
    GenericFunctionUtility& rFunction
    )
{
    rValue = rFunction.RotateAndCallFunction(pEntity->X(), pEntity->Y(), pEntity->Z(), Time, pEntity->X0(), pEntity->Y0(), pEntity->Z0());
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Condition>::CallFunctionLocalSystemComponents(
    const typename Condition::Pointer pEntity,
    const double Time,
    double& rValue,
    GenericFunctionUtility& rFunction
    )
{
    auto& r_entity_geometry = pEntity->GetGeometry();
    const array_1d<double,3>& r_center = r_entity_geometry.Center();
    rValue = rFunction.RotateAndCallFunction(r_center[0],r_center[1],r_center[2], Time);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Element>::CallFunctionLocalSystemComponents(
    const typename Element::Pointer pEntity,
    const double Time,
    double& rValue,
    GenericFunctionUtility& rFunction
    )
{
    auto& r_entity_geometry = pEntity->GetGeometry();
    const array_1d<double,3>& r_center = r_entity_geometry.Center();
    rValue = rFunction.RotateAndCallFunction(r_center[0],r_center[1],r_center[2], Time);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::AssignTimeDependentValue(
    const typename Node::Pointer pEntity,
    const double Time,
    Vector& rValue,
    const double Value
    )
{
    if(rValue.size() !=  1) {
        rValue.resize(1,false);
    }

    rValue[0] = Value;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsHistoricalVariable>::AssignTimeDependentValue(
    const typename Node::Pointer pEntity,
    const double Time,
    Vector& rValue,
    const double Value
    )
{
    if(rValue.size() !=  1) {
        rValue.resize(1,false);
    }

    rValue[0] = Value;
}


/***********************************************************************************/
/***********************************************************************************/

template<>
void AssignScalarFieldToEntitiesProcess<Condition>::AssignTimeDependentValue(
    const typename Condition::Pointer pEntity,
    const double Time,
    Vector& rValue,
    const double Value
    )
{
    auto& r_entity_geometry = pEntity->GetGeometry();
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
    const typename Element::Pointer pEntity,
    const double Time,
    Vector& rValue,
    const double Value
    )
{
    auto& r_entity_geometry = pEntity->GetGeometry();
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

template<class TEntity, bool THistorical>
void AssignScalarFieldToEntitiesProcess<TEntity, THistorical>::InternalAssignValueVector(
    const Variable<Vector>& rVar,
    const double Time
    )
{
    auto& r_entities_array = GetEntitiesContainer();
    const SizeType number_of_entities = r_entities_array.size();

    Vector value;

    if(number_of_entities != 0) {
        auto it_begin = r_entities_array.begin();

        if(mpFunction->DependsOnSpace()) {
            if(mpFunction->UseLocalSystem()) {
                // TODO: Do parallelize eventually (failing for some reason)
                for(IndexType i = 0; i<number_of_entities; ++i) {
                    auto it_entity = it_begin + i;
                    this->CallFunctionLocalSystem(*(it_entity.base()), Time, value, *mpFunction);
                    if constexpr(THistorical == AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable) {
                        it_entity->SetValue(rVar, value);
                    } else {
                        it_entity->FastGetSolutionStepValue(rVar) = value;
                    }
                }
            } else {
                // TODO: Do parallelize eventually (failing for some reason)
                for(IndexType i = 0; i<number_of_entities; ++i) {
                    auto it_entity = it_begin + i;
                    this->CallFunction(*(it_entity.base()), Time, value, *mpFunction);
                    if constexpr(THistorical == AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable) {
                        it_entity->SetValue(rVar, value);
                    } else {
                        it_entity->FastGetSolutionStepValue(rVar) = value;
                    }
                }
            }
        } else { // Only varies in time
            const double time_value = mpFunction->CallFunction(0.0, 0.0, 0.0, Time);
            IndexPartition<IndexType>(number_of_entities).for_each([&](IndexType i) {
                auto it_entity = it_begin + i;
                this->AssignTimeDependentValue(*(it_entity.base()), Time, value, time_value);
                if constexpr(THistorical == AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable) {
                    it_entity->SetValue(rVar, value);
                } else {
                    it_entity->FastGetSolutionStepValue(rVar) = value;
                }
            });
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity, bool THistorical>
void AssignScalarFieldToEntitiesProcess<TEntity, THistorical>::InternalAssignValueScalar(
    const Variable<double>& rVar,
    const double Time
    )
{
    auto& r_entities_array = GetEntitiesContainer();
    const SizeType number_of_entities = r_entities_array.size();

    if(number_of_entities != 0) {
        auto it_begin = r_entities_array.begin();

        if(mpFunction->DependsOnSpace()) {
            double value = 0.0;

            if(mpFunction->UseLocalSystem()) {
                // TODO: Do parallelize eventually (failing for some reason)
                for(IndexType i = 0; i<number_of_entities; ++i) {
                    auto it_entity = it_begin + i;
                    this->CallFunctionLocalSystemComponents(*(it_entity.base()), Time, value, *mpFunction);
                    if constexpr(THistorical == AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable) {
                        it_entity->SetValue(rVar, value);
                    } else {
                        it_entity->FastGetSolutionStepValue(rVar) = value;
                    }
                }
            } else {
                // TODO: Do parallelize eventually (failing for some reason)
                for(IndexType i = 0; i<number_of_entities; ++i) {
                    auto it_entity = it_begin + i;
                    this->CallFunctionComponents(*(it_entity.base()), Time, value, *mpFunction);
                    if constexpr(THistorical == AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable) {
                        it_entity->SetValue(rVar, value);
                    } else {
                        it_entity->FastGetSolutionStepValue(rVar) = value;
                    }
                }
            }
        } else { // Only varies in time
            const double time_value = mpFunction->CallFunction(0.0, 0.0, 0.0, Time);
            IndexPartition<IndexType>(number_of_entities).for_each([&](IndexType i) {
                auto it_entity = it_begin + i;
                if constexpr(THistorical == AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable) {
                    it_entity->SetValue(rVar, time_value);
                } else {
                    it_entity->FastGetSolutionStepValue(rVar) = time_value;
                }
            });
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Node, IndexedObject>& AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable>::GetEntitiesContainer()
{
    return mrModelPart.Nodes();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Node, IndexedObject>& AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsHistoricalVariable>::GetEntitiesContainer()
{
    return mrModelPart.Nodes();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& AssignScalarFieldToEntitiesProcess<Condition>::GetEntitiesContainer()
{
    return mrModelPart.Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>& AssignScalarFieldToEntitiesProcess<Element>::GetEntitiesContainer()
{
    return mrModelPart.Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template class AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsNonHistoricalVariable>;
template class AssignScalarFieldToEntitiesProcess<Node, AssignScalarFieldToEntitiesProcessSettings::SaveAsHistoricalVariable>;
template class AssignScalarFieldToEntitiesProcess<Condition>;
template class AssignScalarFieldToEntitiesProcess<Element>;

}  // namespace Kratos.