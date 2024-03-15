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

#pragma once

// System includes
#include <vector>
#include <type_traits>

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "expression/container_expression.h"
#include "expression/literal_flat_expression.h"
#include "spatial_containers/spatial_containers.h"

// Application includes
#include "entity_point.h"
#include "filter_function.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) ExplicitDampingUtils {
public:
    ///@name Public type definitions
    ///@{

    using IndexType = std::size_t;

    template<class TContainerType>
    using EntityType = typename TContainerType::value_type;

    template<class TContainerType>
    using EntityPointType = EntityPoint<EntityType<TContainerType>>;

    template<class TContainerType>
    using EntityPointsVector = std::vector<typename EntityPointType<TContainerType>::Pointer>;

    // Type definitions for tree-search
    template<class TContainerType>
    using BucketType = Bucket<3, EntityPointType<TContainerType>, EntityPointsVector<TContainerType>>;

    template<class TContainerType>
    using KDTree = Tree<KDTreePartition<BucketType<TContainerType>>>;

    ///@}
    ///@name Public static methods
    ///@{

    template<class TContainerType>
    static ContainerExpression<TContainerType> ComputeDampingCoefficientsBasedOnNearestEntity(
        const ContainerExpression<TContainerType>& rDampingRadiusExpression,
        const std::vector<std::vector<const ModelPart*>>& rDampedModelParts,
        const std::vector<IndexType>& rShape,
        const std::string& rFilterFunctionType,
        const IndexType BucketSize);

    ///@}

private:
    ///@name Private static methods
    ///@{

    template<class TContainerType>
    static void ComputeDampingCoefficientsBasedOnNearestEntityForComponent(
        LiteralFlatExpression<double>& rOutputExpression,
        const ContainerExpression<TContainerType>& rDampingRadiusExpression,
        const std::vector<const ModelPart*>& rDampedModelParts,
        const FilterFunction& rFilterFunction,
        const IndexType BucketSize,
        const IndexType DampedComponentIndex);

    ///@}
};
} // namespace Kratos