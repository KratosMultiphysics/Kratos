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

    using EntityPointVector = std::vector<typename EntityType::Pointer>;

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

    /**
     * @brief Set the damping radius
     */
    virtual void SetRadius(const ContainerExpression<TContainerType>& rExpression)
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::SetRadius. This should be implemented in the derived class.";
    }

    /**
     * @brief Get the used damping radius expression
     */
    virtual typename ContainerExpression<TContainerType>::Pointer GetRadius() const
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::GetRadius. This should be implemented in the derived class.";
        return nullptr;
    }

    /**
     * @brief Get the stride (number of components) this damping can work on
     */
    virtual IndexType GetStride() const
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::GetStride. This should be implemented in the derived class.";
        return 0;
    }

    /**
     * @brief Get the damped model parts
     * @details This method returns list of model parts which are damped for each
     *          component in the stride.
     */
    virtual std::vector<std::vector<ModelPart*>> GetDampedModelParts() const
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::GetDampedModelParts. This should be implemented in the derived class.";
        return {};
    }

    /**
     * @brief Update the damping
     *
     */
    virtual void Update()
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::Update. This should be implemented in the derived class.";
    }

    /**
     * @brief Apply damping to the given weights.
     * @details This method applies damping to the passed @ref rWeights. The damped
     *          weights are stored in the @ref rDampedWeights.
     *
     * @param rDampedWeights            Output damped weights.
     * @param rWeights                  Input filtering weights.
     * @param Index                     Index of the current damping computation.
     * @param NumberOfNeighbours        number of neighbours around the Index point.
     * @param rNeighbours               Neighbour entities around the Index point.
     */
    virtual void Apply(
        std::vector<std::vector<double>>& rDampedWeights,
        const std::vector<double>& rWeights,
        const IndexType Index,
        const IndexType NumberOfNeighbours,
        const EntityPointVector& rNeighbours) const
    {
        KRATOS_ERROR << "Calling base class ExplicitDamping::GetFilteredDampedValue. This should be implemented in the derived class.";
    }

    /**
     * @brief Calculate the damping matrix for the component in stride.
     * @warning Only to be used in testing.
     *
     */
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

