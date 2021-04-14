//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes

// External includes

// Project includes
#include "utilities/sub_model_part_entities_boolean_operation_utility.h"

namespace Kratos
{

template<>
ModelPart::NodesContainerType& SubModelPartEntitiesBooleanOperationUtility<Node<3>,ModelPart::NodesContainerType>::GetContainer(
    ModelPart& rModelPart)
{
    return rModelPart.Nodes();
}

template<>
ModelPart::ElementsContainerType& SubModelPartEntitiesBooleanOperationUtility<Element,ModelPart::ElementsContainerType>::GetContainer(
    ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template<>
ModelPart::ConditionsContainerType& SubModelPartEntitiesBooleanOperationUtility<Condition,ModelPart::ConditionsContainerType>::GetContainer(
    ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

template<>
void SubModelPartEntitiesBooleanOperationUtility<Node<3>,ModelPart::NodesContainerType>::AddEntities(
    const std::vector<IndexType>& rIds,
    ModelPart& rModelPart)
{
    rModelPart.AddNodes(rIds);
}

template<>
void SubModelPartEntitiesBooleanOperationUtility<Element,ModelPart::ElementsContainerType>::AddEntities(
    const std::vector<IndexType>& rIds,
    ModelPart& rModelPart)
{
    rModelPart.AddElements(rIds);
}

template<>
void SubModelPartEntitiesBooleanOperationUtility<Condition,ModelPart::ConditionsContainerType>::AddEntities(
    const std::vector<IndexType>& rIds,
    ModelPart& rModelPart)
{
    rModelPart.AddConditions(rIds);
}

template class SubModelPartEntitiesBooleanOperationUtility<Node<3>,ModelPart::NodesContainerType>;
template class SubModelPartEntitiesBooleanOperationUtility<Element,ModelPart::ElementsContainerType>;
template class SubModelPartEntitiesBooleanOperationUtility<Condition,ModelPart::ConditionsContainerType>;
}  // namespace Kratos.
