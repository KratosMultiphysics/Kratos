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
#include "includes/define.h"
#include "expression/expression.h"
#include "expression/container_expression.h"

// Application includes
#include "entity_point.h"

namespace Kratos
{

///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) ExplicitDamping
{
  public:
    ///@name Type definitions
    ///@{

    using EntityType = typename TContainerType::value_type;

    using EntityPointType = EntityPoint<EntityType>;

    using EntityPointVector = std::vector<typename EntityPointType::Pointer>;

    /// Pointer definition of ContainerMapper
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitDamping);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ExplicitDamping() = default;

    virtual ~ExplicitDamping() = default;

    ///@}
    ///@name Operations
    ///@{

    virtual void SetRadius(const ContainerExpression<TContainerType>& rExpression)
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::SetRadius. This should be implemented in the derived class.";
    }

    virtual typename ContainerExpression<TContainerType>::Pointer GetRadius() const
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::GetRadius. This should be implemented in the derived class.";
        return nullptr;
    }

    virtual IndexType GetStride() const
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::GetStride. This should be implemented in the derived class.";
        return 0;
    }

    virtual std::vector<std::vector<ModelPart*>> GetDampedModelParts() const
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::GetDampedModelParts. This should be implemented in the derived class.";
        return {};
    }

    virtual void Update()
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::Update. This should be implemented in the derived class.";
    }

    virtual void Apply(
        std::vector<std::vector<double>>& rDampedWeights,
        const std::vector<double>& rWeights,
        const IndexType Index,
        const IndexType NumberOfNeighbours,
        const EntityPointVector& rNeighbours) const
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::GetFilteredDampedValue. This should be implemented in the derived class.";
    }

    virtual void CalculateMatrix(
        Matrix& rOutput,
        const IndexType ComponentIndex) const
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::CalculateMatrix. This should be implemented in the derived class.";
    }

    ///@}

}; // Class ExplicitDamping

///@}

} // namespace Kratos.

