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

// Project includes
#include "containers/model.h"
#include "containers/pointer_vector_set.h"
#include "expression/container_expression.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"

// Application includes
#include "entity_point.h"
#include "filter.h"
#include "filter_function.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) ExplicitFilter: public Filter<TContainerType>
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
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitFilter);

    ///@}
    ///@name LifeCycle
    ///@{

    ExplicitFilter(
        Model& rModel,
        Parameters Settings);

    ~ExplicitFilter() override = default;

    ///@}
    ///@name Public operations

    void Initialize() override;

    void Check() override {};

    void SetFilterRadius(const ContainerExpression<TContainerType>& rContainerExpression);

    void Update() override;

    void Finalize() override {}

    ContainerExpression<TContainerType> FilterField(const ContainerExpression<TContainerType>& rContainerExpression) const override;

    ContainerExpression<TContainerType> FilterIntegratedField(const ContainerExpression<TContainerType>& rContainerExpression) const override;

    void GetIntegrationWeights(ContainerExpression<TContainerType>& rContainerExpression) const;

    std::string Info() const override;

    ///@}
private:
    ///@name Private member variables
    ///@{

    Model& mrModel;

    std::string mFilteringModelPartName;

    FilterFunction::UniquePointer mpKernelFunction;

    FilterFunction::UniquePointer mpDampingFunction;

    typename ContainerExpression<TContainerType>::Pointer mpFilterRadiusContainer;

    Expression::ConstPointer mpNodalDomainSizeExpression;

    EntityPointVector mEntityPointVector;

    IndexType mBucketSize = 100;

    IndexType mMaxNumberOfNeighbors;

    IndexType mNumberOfComponents;

    double mFilterRadius;

    std::vector<std::string> mDampingModelPartNames;

    std::vector<std::vector<bool>> mDampingComponentIndices;

    typename KDTree::Pointer mpSearchTree;

    ///@}
    ///@name Private operations
    ///@{

    template<class TWeightIntegrationType>
    ContainerExpression<TContainerType> GenericFilterField(const ContainerExpression<TContainerType>& rContainerExpression) const;

    ///@}
};

///@}
} // namespace Kratos