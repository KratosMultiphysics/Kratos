//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//                   Suneth Warnakulasuriya
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
#include "damping_function.h"
#include "explicit_damping.h"
#include "filter_utils.h"

namespace Kratos
{

///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) IntegratedNearestEntityExplicitDamping: public ExplicitDamping<TContainerType>
{
  public:
    ///@name Type definitions
    ///@{

    using BaseType = ExplicitDamping<TContainerType>;

    using EntityType = typename BaseType::EntityType;

    using EntityPointVector = typename BaseType::EntityPointVector;

    using PositionAdapter = NanoFlannMultipleModelPartPositionAdapter<TContainerType>;

    using ResultVectorType = typename PositionAdapter::ResultVectorType;

    using PointerVectorType = typename PositionAdapter::PointerVectorType;

    using DistanceMetricType = typename nanoflann::metric_L2_Simple::traits<double, PositionAdapter>::distance_t;

    using KDTreeIndexType = nanoflann::KDTreeSingleIndexAdaptor<DistanceMetricType, PositionAdapter, 3>;

    using KDTreeThreadLocalStorage = NanoFlannKDTreeThreadLocalStorage<PointerVectorType>;

    /// Pointer definition of ContainerMapper
    KRATOS_CLASS_POINTER_DEFINITION(IntegratedNearestEntityExplicitDamping);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    IntegratedNearestEntityExplicitDamping(
        Model& rModel,
        Parameters Settings,
        const IndexType Stride);

    virtual ~IntegratedNearestEntityExplicitDamping() = default;

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

    IndexType mLeafMaxSize;

    DampingFunction::UniquePointer mpKernelFunction;

    typename ContainerExpression<TContainerType>::Pointer mpDampingRadius;

    std::vector<std::vector<ModelPart*>> mComponentWiseDampedModelParts;

    std::vector<std::shared_ptr<PositionAdapter>> mComponentWisePositionAdapters;

    std::vector<std::shared_ptr<KDTreeIndexType>> mComponentWiseKDTrees;

    ///@}

}; // Class IntegratedNearestEntityExplicitDamping

///@}

} // namespace Kratos.

