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
#include "utilities/parallel_utilities.h"
#include "processes/entity_erase_process.h"
#include "utilities/sub_model_part_entities_boolean_operation_utility.h"

namespace Kratos
{

template<>
ModelPart::NodesContainerType& SubModelPartEntitiesBooleanOperationUtility<
    Node<3>,ModelPart::NodesContainerType>::GetContainer(ModelPart& rModelPart)
{
    return rModelPart.Nodes();
}

template<>
ModelPart::ElementsContainerType& SubModelPartEntitiesBooleanOperationUtility<
    Element,ModelPart::ElementsContainerType>::GetContainer(ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template<>
ModelPart::ConditionsContainerType& SubModelPartEntitiesBooleanOperationUtility<
    Condition,ModelPart::ConditionsContainerType>::GetContainer(ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

template<>
ModelPart::MasterSlaveConstraintContainerType& SubModelPartEntitiesBooleanOperationUtility<
    MasterSlaveConstraint,ModelPart::MasterSlaveConstraintContainerType>::GetContainer(ModelPart& rModelPart)
{
    return rModelPart.MasterSlaveConstraints();
}

template<>
void SubModelPartEntitiesBooleanOperationUtility<
    Node<3>,ModelPart::NodesContainerType>::AddEntities(
        const std::vector<IndexType>& rIds,ModelPart& rModelPart)
{
    rModelPart.AddNodes(rIds);
}

template<>
void SubModelPartEntitiesBooleanOperationUtility<
    Element,ModelPart::ElementsContainerType>::AddEntities(
        const std::vector<IndexType>& rIds,
        ModelPart& rModelPart)
{
    rModelPart.AddElements(rIds);
}

template<>
void SubModelPartEntitiesBooleanOperationUtility<
    Condition,ModelPart::ConditionsContainerType>::AddEntities(
        const std::vector<IndexType>& rIds,
        ModelPart& rModelPart)
{
    rModelPart.AddConditions(rIds);
}

template<>
void SubModelPartEntitiesBooleanOperationUtility<
    MasterSlaveConstraint,ModelPart::MasterSlaveConstraintContainerType>::AddEntities(
        const std::vector<IndexType>& rIds,
        ModelPart& rModelPart)
{
    rModelPart.AddMasterSlaveConstraints(rIds);
}

template<class TEntityType, class TContainerType>
std::vector<IndexType> SubModelPartEntitiesBooleanOperationUtility<
    TEntityType, TContainerType>::GetContainerIds(ModelPart& rModelPart)
{
    const TContainerType& r_container = GetContainer(rModelPart);
    std::vector<IndexType> ids_vector(r_container.size());
    IndexPartition<std::size_t>(r_container.size()).for_each([&](std::size_t i){
        ids_vector[i] = (r_container.begin()+i)->Id();
    });
    return ids_vector;
}

template<class TEntityType, class TContainerType>
void SubModelPartEntitiesBooleanOperationUtility<
    TEntityType, TContainerType>::BooleanOperation(
    ModelPart& rModelPartA,
    ModelPart& rModelPartB,
    ModelPart& rDestination,
    BooleanOperators ThisOperator)
{
    KRATOS_ERROR_IF(!rDestination.IsSubModelPart()) << "The destination model part must be a sub model part." << std::endl;
    const ModelPart& r_root_a = rModelPartA.GetRootModelPart();
    const ModelPart& r_root_b = rModelPartB.GetRootModelPart();
    const ModelPart& r_root_d = rDestination.GetRootModelPart();
    KRATOS_ERROR_IF(&r_root_a != &r_root_b) << "The first and second model parts must belong to the same root model part." << std::endl;
    KRATOS_ERROR_IF(&r_root_a != &r_root_d) << "The destination model part must belong to the same root model part than the first and the second." << std::endl;
    std::vector<IndexType> ids_a = GetContainerIds(rModelPartA);
    std::vector<IndexType> ids_b = GetContainerIds(rModelPartB);
    std::vector<IndexType> ids_destination;
    std::sort(ids_a.begin(), ids_a.end());
    std::sort(ids_b.begin(), ids_b.end());

    if (ThisOperator == BooleanOperators::Union) {
        std::set_union(
            ids_a.begin(), ids_a.end(),
            ids_b.begin(), ids_b.end(),
            std::back_inserter(ids_destination));
    } else if (ThisOperator == BooleanOperators::Intersection) {
        std::set_intersection(
            ids_a.begin(), ids_a.end(),
            ids_b.begin(), ids_b.end(),
            std::back_inserter(ids_destination));
    } else if (ThisOperator == BooleanOperators::Difference) {
        std::set_difference(
            ids_a.begin(), ids_a.end(),
            ids_b.begin(), ids_b.end(),
            std::back_inserter(ids_destination));
    }

    EntitiesEraseProcess<TEntityType>(rDestination, EntitiesEraseProcessFlags::ERASE_ALL_ENTITIES).Execute();
    AddEntities(ids_destination, rDestination);
}

template class SubModelPartEntitiesBooleanOperationUtility<Node<3>,ModelPart::NodesContainerType>;
template class SubModelPartEntitiesBooleanOperationUtility<Element,ModelPart::ElementsContainerType>;
template class SubModelPartEntitiesBooleanOperationUtility<Condition,ModelPart::ConditionsContainerType>;
template class SubModelPartEntitiesBooleanOperationUtility<MasterSlaveConstraint,ModelPart::MasterSlaveConstraintContainerType>;
}  // namespace Kratos.
