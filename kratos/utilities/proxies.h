//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// Project includes
#include "includes/global_variables.h" // DataLocation
#include "includes/node.h" // Node
#include "includes/element.h" // Element
#include "includes/condition.h" // Condition

// System includes
#include <type_traits> // remove_reference_t, is_const_v, is_same_v, decay_t
#include <optional> // optional


namespace Kratos {


template <Globals::DataLocation TLocation, bool TMutable>
class EntityProxy
{
private:
    constexpr static Globals::DataLocation Location = TLocation;

    /// @ref Node, @ref Element, or @ref Condition without a const-qualifier, depending on @a TLocation.
    using UnqualifiedEntity = std::conditional_t<
        TLocation == Globals::DataLocation::NodeHistorical || TLocation == Globals::DataLocation::NodeNonHistorical,
        Node,
        std::conditional_t<
            TLocation == Globals::DataLocation::Element,
            Element,
            std::conditional_t<
                TLocation == Globals::DataLocation::Condition,
                Condition,
                void // <== invalid fallback type; will throw a compile-time error
            >
        >
    >;

    /// Const-qualified @ref Node, @ref Element, or @ref Condition, depending on @a TLocation and @a TMutable.
    using QualifiedEntity = std::conditional_t<TMutable,
                                               UnqualifiedEntity,
                                               const UnqualifiedEntity>;

    /// @ref ContainerProxy needs to access private typedefs.
    friend class ContainerProxy<EntityProxy>;

public:
    /// @brief Default constructor that leaves the instance in an invalid state.
    /// @throws if any member function is called without reassigning this instance
    ///         with a valid one.
    EntityProxy() noexcept = default;

    /// @brief Constructor creating a valid proxy, wrapping the input entity.
    /// @param rEntity Entity that will be accessed when member functions are called.
    /// @warning This proxy is invalidated when the container holding @a rEntity
    ///          invalidates its iterators or when @a rEntity is destroyed.
    EntityProxy(QualifiedEntity& rEntity) noexcept : mpEntity(&rEntity) {}

    /// @brief Fetch the value corresponding to the input variable in the wrapped entity.
    template <class TVariable>
    const typename TVariable::Type& GetValue(const TVariable& rVariable) const
    {
        if constexpr (TLocation == Globals::DataLocation::NodeHistorical) {
            return mpEntity.value().GetValue(rVariable);
        } else {
            return mpEntity.value().GetSolutionStepValue(rVariable);
        }
    }

    /// @brief Fetch the value corresponding to the input variable in the wrapped entity.
    template <class TVariable, std::enable_if_t<TMutable> = true>
    typename TVariable::Type& GetValue(const TVariable& rVariable)
    {
        if constexpr (TLocation == Globals::DataLocation::NodeHistorical) {
            return mpEntity.value().GetValue(rVariable);
        } else {
            return mpEntity.value().GetSolutionStepValue(rVariable);
        }
    }

    /// @brief Overwrite the value corresponding to the input variable in the wrapped entity.
    template <class TVariable, std::enable_if_t<TMutable,bool> = true>
    void SetValue(const TVariable& rVariable, const typename TVariable::Type& rValue)
    {
        mpEntity.value().SetValue(rVariable, rValue);
    }

    /// @brief Immutable access to the wrapped entity.
    const UnqualifiedEntity& GetEntity() const
    {
        return *mpEntity.value();
    }

    /// @brief Mutable or immutable access to the wrapped entity, depending on @a TMutable.
    QualifiedEntity& GetEntity()
    {
        return *mpEntity.value();
    }

private:
    std::optional<QualifiedEntity*> mpEntity;
}; // class EntityProxy


} // namespace Kratos
