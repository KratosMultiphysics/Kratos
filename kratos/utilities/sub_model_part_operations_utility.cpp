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
#include "utilities/sub_model_part_operations_utility.h"

namespace Kratos
{

template<>
ModelPart::NodesContainerType& SubModelPartOperationsUtility::GetContainer<ModelPart::NodesContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Nodes();
}

template<>
ModelPart::ElementsContainerType& SubModelPartOperationsUtility::GetContainer<ModelPart::ElementsContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template<>
ModelPart::ConditionsContainerType& SubModelPartOperationsUtility::GetContainer<ModelPart::ConditionsContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

template<>
void SubModelPartOperationsUtility::AddEntities<Node<3>>(
    const std::vector<IndexType>& rIds,
    ModelPart& rModelPart)
{
    rModelPart.AddNodes(rIds);
}

template<>
void SubModelPartOperationsUtility::AddEntities<Element>(
    const std::vector<IndexType>& rIds,
    ModelPart& rModelPart)
{
    rModelPart.AddElements(rIds);
}

template<>
void SubModelPartOperationsUtility::AddEntities<Condition>(
    const std::vector<IndexType>& rIds,
    ModelPart& rModelPart)
{
    rModelPart.AddConditions(rIds);
}

}  // namespace Kratos.
