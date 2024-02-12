//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <string>
#include <type_traits>

// Project includes
#include "includes/define.h"
#include "expression/container_expression.h"
#include "expression/domain_size_expression_io.h"

// Application includes
#include "custom_utilities/container_expression_utils.h"

// Include base h
#include "entity_node_entity_filter.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TContainerType>
EntityNodeEntityFilter<TContainerType>::EntityNodeEntityFilter(ModelPart& rModelPart)
    : mrModelPart(rModelPart),
      mNeighbourEntities(rModelPart),
      mEntityDomainSize(rModelPart)
{
    Update();
}

template<class TContainerType>
void EntityNodeEntityFilter<TContainerType>::Update()
{
    ContainerExpressionUtils::ComputeNumberOfNeighbourEntities<TContainerType>(mNeighbourEntities);
    DomainSizeExpressionIO::Read(mEntityDomainSize);
}

template<class TContainerType>
ContainerExpression<TContainerType> EntityNodeEntityFilter<TContainerType>::FilterField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    // for conditions and elements, first we map the field to nodes, and map them back to respective conditions or elements
    if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType> || std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        ContainerExpression<ModelPart::NodesContainerType> nodal_expression(*rContainerExpression.pGetModelPart());
        ContainerExpressionUtils::MapContainerVariableToNodalVariable(nodal_expression, rContainerExpression, mNeighbourEntities);

        ContainerExpression<TContainerType> output(*rContainerExpression.pGetModelPart());
        ContainerExpressionUtils::MapNodalVariableToContainerVariable(output, nodal_expression);
        return output;
    } else {
        static_assert(std::is_same_v<TContainerType, TContainerType>, "Unsupported TContainerType");
    }
}

template<class TContainerType>
ContainerExpression<TContainerType> EntityNodeEntityFilter<TContainerType>::FilterIntegratedField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    // first divide by the entity area
    if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType> || std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        auto copy = rContainerExpression;
        return FilterField(copy / mEntityDomainSize);
    } else {
        static_assert(std::is_same_v<TContainerType, TContainerType>, "Unsupported TContainerType");
    }
}

template<class TContainerType>
std::string EntityNodeEntityFilter<TContainerType>::Info() const
{
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        return "NodeEntityNodeFilter";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        return "ConditionNodeConditionFilter";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        return "ElementNodeElementFilter";
    } else {
        static_assert(std::is_same_v<TContainerType, TContainerType>, "Unsupported TContainerType");
    }
}

// template instantiations
// template class EntityNodeEntityFilter<ModelPart::ConditionsContainerType>;
template class EntityNodeEntityFilter<ModelPart::ConditionsContainerType>;
template class EntityNodeEntityFilter<ModelPart::ElementsContainerType>;
} // namespace Kratos