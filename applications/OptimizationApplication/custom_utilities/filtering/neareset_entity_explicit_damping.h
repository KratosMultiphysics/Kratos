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

// Project includes
#include "containers/model.h"
#include "expression/container_expression.h"
#include "expression/expression.h"
#include "includes/define.h"
#include "spatial_containers/spatial_containers.h"

// Application includes
#include "filter_function.h"
#include "explicit_damping.h"

namespace Kratos
{

///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) NearestEntityExplicitDamping: public ExplicitDamping<TContainerType>
{
  public:
    ///@name Type definitions
    ///@{

    using BaseType = ExplicitDamping<TContainerType>;

    using EntityType = typename BaseType::EntityType;

    using EntityPointType = typename BaseType::EntityPointType;

    using EntityPointVector = typename BaseType::EntityPointVector;

    // Type definitions for tree-search
    using BucketType = Bucket<3, EntityPointType, EntityPointVector>;

    using KDTree = Tree<KDTreePartition<BucketType>>;

    /// Pointer definition of ContainerMapper
    KRATOS_CLASS_POINTER_DEFINITION(NearestEntityExplicitDamping);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    NearestEntityExplicitDamping(
        Model& rModel,
        Parameters Settings,
        const IndexType Stride);

    virtual ~NearestEntityExplicitDamping() = default;

    ///@}
    ///@name Operations
    ///@{

    void SetRadius(const ContainerExpression<TContainerType>& rDampingRadiusExpression) override;

    typename ContainerExpression<TContainerType>::Pointer GetRadius() const override;

    IndexType GetStride() const override;

    std::vector<std::vector<ModelPart*>> GetDampedModelParts() const override;

    virtual void Apply(
        std::vector<std::vector<double>>& rDampedWeights,
        const std::vector<double>& rWeights,
        const IndexType Index,
        const IndexType NumberOfNeighbours,
        const EntityPointVector& rNeighbours) const override;

    void Update() override;

    void CalculateMatrix(
        Matrix& rOutput,
        const IndexType ComponentIndex) const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    IndexType mStride;

    IndexType mBucketSize = 100;

    FilterFunction::UniquePointer mpKernelFunction;

    Expression::ConstPointer mpDampingCoefficients;

    typename ContainerExpression<TContainerType>::Pointer mpDampingRadius;

    std::vector<std::vector<ModelPart*>> mComponentWiseDampedModelParts;

    ///@}

}; // Class NearestEntityExplicitDamping

///@}

} // namespace Kratos.

