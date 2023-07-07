//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel
//                   Aditya Ghantasala
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>
#include <variant>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "expression/container_expression.h"

// Application includes
#include "entity_point.h"
#include "filter_function.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) ExplicitVertexMorphingFilter
{
public:
    ///@name Type definitions
    ///@{

    using EntityType = typename TContainerType::value_type;

    using EntityPointVector = std::vector<typename EntityPoint<EntityType>::Pointer>;

    // Type definitions for tree-search
    using BucketType = Bucket<3, EntityPoint<EntityType>, EntityPointVector>;

    using KDTree = Tree<KDTreePartition<BucketType>>;

    /// Pointer definition of ContainerMapper
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitVertexMorphingFilter);

    ///@}
    ///@name LifeCycle
    ///@{

    ExplicitVertexMorphingFilter(
        const ModelPart& rModelPart,
        const std::string& rKernelFunctionType,
        const IndexType MaxNumberOfNeighbours);

    ///@}
    ///@name Public operations

    void SetFilterRadius(const ContainerExpression<TContainerType>& rContainerExpression);

    void Update();

    ContainerExpression<TContainerType> FilterField(const ContainerExpression<TContainerType>& rContainerExpression) const;

    ContainerExpression<TContainerType> FilterIntegratedField(const ContainerExpression<TContainerType>& rContainerExpression) const;

    std::string Info() const;

    ///@}
private:
    ///@name Private member variables
    ///@{

    const ModelPart& mrModelPart;

    FilterFunction::UniquePointer mpKernelFunction;

    typename ContainerExpression<TContainerType>::Pointer mpFilterRadiusContainer;

    Expression::ConstPointer mpNodalDomainSizeExpression;

    EntityPointVector mEntityPointVector;

    IndexType mBucketSize = 100;

    IndexType mMaxNumberOfNeighbors;

    typename KDTree::Pointer mpSearchTree;

    ///@}
    ///@name Private operations
    ///@{

    template<class TWeightIntegrationType>
    ContainerExpression<TContainerType> GenericFilterField(const ContainerExpression<TContainerType>& rContainerExpression) const;

    ///@}
};

///@name Input and output
///@{

/// output stream function
template<class TContainerType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ExplicitVertexMorphingFilter<TContainerType>& rThis)
{
    return rOStream << rThis.Info();
}

///@}

///@}
} // namespace Kratos