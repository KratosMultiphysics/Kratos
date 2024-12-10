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
#include "includes/model_part.h" // ModelPart::NodesContainerType, ModelPart::ElementsContainerType, ModelPart::ConditionsContainerType
#include "utilities/variable_utils.h" // VariableUtils::HasValue, VariableUtils::GetValue, VariableUtils::SetValue
#include "utilities/model_part_utils.h" // ModelPartUtils::GetContainer

// System includes
#include <type_traits> // remove_reference_t, is_const_v, is_same_v, decay_t


namespace Kratos {


template <class TEntityProxy>
class ContainerProxy;



/** @brief Wrapper class providing a uniform interface for historical/non-historical @ref Node, @ref Element, and @ref Condition.
 *  @details @ref EntityProxy exposes common functionality related accessing stored @ref Variable s within an entity,
 *         without additional runtime overhead. In this context, an entity can refer to:
 *         - a @ref Node with historical variables (@ref Globals::DataLocation::NodeHistorical)
 *         - a @ref Node with non-historical variables (@ref Globals::DataLocation::NodeNonHistorical)
 *         - an @ref Element (@ref Globals::DataLocation::Element)
 *         - a @ref Condition (@ref Globals::DataLocation::Condition)
 *         The exposed common functionalities include checking, reading and overwriting the values
 *         related to the provided variables associated with the entity.
 *  @warning Default constructed @ref EntityProxy instances are in an invalid state,
 *           and their member functions must not be called.
 *  @throws if member functions of a default constructed instance are called.
 */
template <Globals::DataLocation TLocation, bool TMutable>
class EntityProxy
{
private:
    constexpr static Globals::DataLocation Location = TLocation;

    constexpr static bool IsMutable = TMutable;

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
                std::conditional_t<
                    TLocation == Globals::DataLocation::Constraint,
                    MasterSlaveConstraint,
                    std::conditional_t<
                        TLocation == Globals::DataLocation::ProcessInfo,
                        ProcessInfo,
                        std::conditional_t<
                            TLocation == Globals::DataLocation::ModelPart,
                            ModelPart,
                            void // <== invalid fallback type; will throw a compile-time error
                        >
                    >
                >
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
    /// @brief Delete the default constructor to prevent storing nullptr.
    EntityProxy() = delete;

    /// @brief Constructor creating a valid proxy, wrapping the input entity.
    /// @param rEntity Entity that will be accessed when member functions are called.
    /// @warning This proxy is invalidated when the container holding @a rEntity
    ///          invalidates its iterators or when @a rEntity is destroyed.
    EntityProxy(QualifiedEntity& rEntity) noexcept : mpEntity(&rEntity) {}

    /// @brief Check whether the entity has a value for the provided variable.
    template <class TValue>
    bool HasValue(const Variable<TValue>& rVariable) const noexcept
    {
        return VariableUtils::HasValue<TLocation>(*mpEntity, rVariable);
    }

    /// @brief Fetch the value corresponding to the input variable in the wrapped entity.
    template <class TValue>
    std::conditional_t<std::is_integral_v<TValue> || std::is_floating_point_v<TValue>,
                       TValue,           // <== return by value if scalar type
                       const TValue&>    // <== return by reference in non-scalar type
    GetValue(const Variable<TValue>& rVariable) const
    {
        return VariableUtils::GetValue<TLocation>(*mpEntity, rVariable);
    }

    /// @brief Fetch the value corresponding to the input variable in the wrapped entity.
    template <class TValue, std::enable_if_t</*this is required for SFINAE*/!std::is_same_v<TValue,void> && TMutable,bool> = true>
    TValue& GetValue(const Variable<TValue>& rVariable)
    {
        return VariableUtils::GetValue<TLocation>(*mpEntity, rVariable);
    }

    /// @brief Overwrite the value corresponding to the input variable in the wrapped entity.
    template <class TValue, std::enable_if_t</*this is required for SFINAE*/!std::is_same_v<TValue,void> && TMutable,bool> = true>
    void SetValue(const Variable<TValue>& rVariable,
                  std::conditional_t<std::is_integral_v<TValue> || std::is_floating_point_v<TValue>,
                                     TValue,         /*pass scalar types by value*/
                                     const TValue&>  /*pass non-scalar types by reference*/ Value)
    {
        VariableUtils::SetValue<TLocation>(*mpEntity, rVariable, Value);
    }

    /// @brief Immutable access to the wrapped entity.
    const UnqualifiedEntity& GetEntity() const
    {
        return *mpEntity;
    }

    /// @brief Mutable or immutable access to the wrapped entity, depending on @a TMutable.
    QualifiedEntity& GetEntity()
    {
        return *mpEntity;
    }

private:
    QualifiedEntity* mpEntity;
}; // class EntityProxy



/** @brief A view with a uniform interface for @ref ModelPart::NodesContainerType, @ref ModelPart::ElementsContainerType, or @ref ModelPart::ConditionsContainerType.
 *  @details @ref ContainerProxy provides uniform access to the @ref Variable s stored in the entities of a container.
 *           Entities in the container are wrapped by @ref EntityProxy. In this context, an entity can refer to:
 *         - a @ref Node with historical variables (@ref Globals::DataLocation::NodeHistorical)
 *         - a @ref Node with non-historical variables (@ref Globals::DataLocation::NodeNonHistorical)
 *         - an @ref Element (@ref Globals::DataLocation::Element)
 *         - a @ref Condition (@ref Globals::DataLocation::Condition)
 */
template <class TEntityProxy>
class ContainerProxy
{
private:
    using UnqualifiedContainer = std::conditional_t<
        std::is_same_v<typename TEntityProxy::UnqualifiedEntity,Node>,
        ModelPart::NodesContainerType,
        std::conditional_t<
            std::is_same_v<typename TEntityProxy::UnqualifiedEntity,Element>,
            ModelPart::ElementsContainerType,
            std::conditional_t<
                std::is_same_v<typename TEntityProxy::UnqualifiedEntity,Condition>,
                ModelPart::ConditionsContainerType,
                std::conditional_t<
                    std::is_same_v<typename TEntityProxy::UnqualifiedEntity,MasterSlaveConstraint>,
                    ModelPart::MasterSlaveConstraintContainerType,
                    void // <== invalid fallback type; will throw a compile-time error
                >
            >
        >
    >;

    constexpr static bool IsMutable = TEntityProxy::IsMutable;

    using WrappedIterator = std::conditional_t<IsMutable,
                                               typename UnqualifiedContainer::iterator,
                                               typename UnqualifiedContainer::const_iterator>;

    template <bool TMutable>
    class Iterator
    {
    private:
        using Wrapped = std::conditional_t<TMutable,
                                           typename UnqualifiedContainer::iterator,
                                           typename UnqualifiedContainer::const_iterator>;

    public:
        using value_type = EntityProxy<TEntityProxy::Location,TMutable>;

        using pointer = std::conditional_t<TMutable,
                                           value_type*,
                                           const value_type*>;

        using reference = std::conditional_t<TMutable,
                                             value_type&,
                                             const value_type&>;

        using difference_type = std::ptrdiff_t;

        using iterator_category = std::random_access_iterator_tag;

        Iterator() noexcept = default;

        Iterator(Wrapped It) noexcept : mWrapped(It) {}

        value_type operator*() const noexcept {return value_type(*mWrapped);}

        Iterator& operator++() noexcept {++mWrapped; return *this;}

        Iterator operator++(int) noexcept {Iterator copy(mWrapped); ++mWrapped; return copy;}

        Iterator& operator--() noexcept {--mWrapped; return *this;}

        Iterator operator--(int) noexcept {Iterator copy(mWrapped); --mWrapped; return copy;}

        Iterator& operator+=(difference_type Rhs) noexcept {mWrapped += Rhs; return *this;}

        Iterator& operator-=(difference_type Rhs) noexcept {mWrapped -= Rhs; return *this;}

        Iterator operator+(difference_type Rhs) const noexcept {Iterator copy(mWrapped); copy += Rhs; return copy;}

        Iterator operator-(difference_type Rhs) const noexcept {Iterator copy(mWrapped); copy -= Rhs; return copy;}

        difference_type operator-(Iterator Rhs) const noexcept {return mWrapped - Rhs.mWrapped;}

        bool operator==(Iterator Rhs) const noexcept {return mWrapped == Rhs.mWrapped;}

        bool operator!=(Iterator Rhs) const noexcept {return mWrapped != Rhs.mWrapped;}

        bool operator<(Iterator Rhs) const noexcept {return mWrapped < Rhs.mWrapped;}

        bool operator>(Iterator Rhs) const noexcept {return mWrapped > Rhs.mWrapped;}

        bool operator<=(Iterator Rhs) const noexcept {return mWrapped <= Rhs.mWrapped;}

        bool operator>=(Iterator Rhs) const noexcept {return mWrapped >= Rhs.mWrapped;}

    private:
        Wrapped mWrapped;
    }; // class Iterator
public:
    using iterator = Iterator<IsMutable>;

    using const_iterator = Iterator<false>;

    using size_type = std::size_t;

    using value_type = typename iterator::value_type;

    ContainerProxy() noexcept = default;

    ContainerProxy(WrappedIterator Begin, WrappedIterator End) noexcept
        : mBegin(Begin),
          mEnd(End)
    {}

    typename const_iterator::value_type operator[](size_type Index) const noexcept {return typename const_iterator::value_type(*(mBegin + Index));}

    typename iterator::value_type operator[](size_type Index) noexcept {return typename iterator::value_type(*(mBegin + Index));}

    typename const_iterator::value_type at(size_type Index) const noexcept {return typename const_iterator::value_type(*(mBegin + Index));}

    typename iterator::value_type at(size_type Index) noexcept {return typename iterator::value_type(*(mBegin + Index));}

    size_type size() const noexcept {return std::distance(mBegin, mEnd);}

    bool empty() const noexcept {return this->size() == 0;}

    const_iterator cbegin() const noexcept {return const_iterator(mBegin);}

    const_iterator begin() const noexcept {return this->cbegin();}

    iterator begin() noexcept {return iterator(mBegin);}

    const_iterator cend() const noexcept {return const_iterator(mEnd);}

    const_iterator end() const noexcept {return this->cend();}

    iterator end() noexcept {return iterator(mEnd);}

private:
    WrappedIterator mBegin, mEnd;
}; // class ContainerProxy



/// @brief Invalid template base to be specialized for valid template parameters.
template <Globals::DataLocation TLocation, class TEntity>
inline auto MakeProxy(const TEntity& rEntity)
{
    static_assert(std::is_same_v<TEntity,void>, "Invalid DataLocation-Entity combination");
}


/// @brief Invalid template base to be specialized for valid template parameters.
template <Globals::DataLocation TLocation, class TEntity>
inline auto MakeProxy(TEntity& rEntity)
{
    static_assert(std::is_same_v<TEntity,void>, "Invalid DataLocation-Entity combination");
}


#define KRATOS_DEFINE_ENTITY_PROXY_FACTORY(TLocation, TEntity)                                  \
    /** @brief Convenience function for constructing immutable @ref EntityProxy instances.*/    \
    template <>                                                                                 \
    inline auto MakeProxy<TLocation,TEntity>(const TEntity& rEntity)    \
    {return EntityProxy<TLocation,false>(rEntity);}                                             \
    /** @brief Convenience function for constructing mutable @ref EntityProxy instances.*/      \
    template <>                                                                                 \
    inline auto MakeProxy<TLocation,TEntity>(TEntity& rEntity)           \
    {return EntityProxy<TLocation,true>(rEntity);}

KRATOS_DEFINE_ENTITY_PROXY_FACTORY(Globals::DataLocation::NodeHistorical, Node)

KRATOS_DEFINE_ENTITY_PROXY_FACTORY(Globals::DataLocation::NodeNonHistorical, Node)

KRATOS_DEFINE_ENTITY_PROXY_FACTORY(Globals::DataLocation::Element, Element)

KRATOS_DEFINE_ENTITY_PROXY_FACTORY(Globals::DataLocation::Condition, Condition)

KRATOS_DEFINE_ENTITY_PROXY_FACTORY(Globals::DataLocation::Constraint, MasterSlaveConstraint)

KRATOS_DEFINE_ENTITY_PROXY_FACTORY(Globals::DataLocation::ProcessInfo, ProcessInfo)

KRATOS_DEFINE_ENTITY_PROXY_FACTORY(Globals::DataLocation::ModelPart, ModelPart)

#undef KRATOS_DEFINE_ENTITY_PROXY_FACTORY


/// @brief Convenience function for constructing a mutable @ref ProcessInfo proxy from a @ref ModelPart.
template <>
inline auto MakeProxy<Globals::DataLocation::ProcessInfo,ModelPart>(const ModelPart& rModelPart)
{
    return EntityProxy<Globals::DataLocation::ProcessInfo,false>(rModelPart.GetProcessInfo());
}


/// @brief Convenience function for constructing an immutable @ref ProcessInfo proxy from a @ref ModelPart.
template <>
inline auto MakeProxy<Globals::DataLocation::ProcessInfo,ModelPart>(ModelPart& rModelPart)
{
    return EntityProxy<Globals::DataLocation::ProcessInfo,true>(rModelPart.GetProcessInfo());
}


#define KRATOS_DEFINE_CONTAINER_PROXY_FACTORY(TLocation)                                        \
    /** @brief Convenience function for constructing immutable @ref ContainerProxy instances.*/ \
    template <>                                                                                 \
    inline auto MakeProxy<TLocation,ModelPart>(const ModelPart& rModelPart)                     \
    {                                                                                           \
        const auto& r_container = ModelPartUtils::GetContainer<TLocation>(rModelPart);          \
        return ContainerProxy<EntityProxy<TLocation,false>>(r_container.begin(),                \
                                                            r_container.end());                 \
    }                                                                                           \
    /** @brief Convenience function for constructing mutable @ref ContainerProxy instances.*/   \
    template <>                                                                                 \
    inline auto MakeProxy<TLocation,ModelPart>(ModelPart& rModelPart)                           \
    {                                                                                           \
        auto& r_container = ModelPartUtils::GetContainer<TLocation>(rModelPart);                \
        return ContainerProxy<EntityProxy<TLocation,true>>(r_container.begin(),                 \
                                                           r_container.end());                  \
    }

KRATOS_DEFINE_CONTAINER_PROXY_FACTORY(Globals::DataLocation::NodeHistorical)

KRATOS_DEFINE_CONTAINER_PROXY_FACTORY(Globals::DataLocation::NodeNonHistorical)

KRATOS_DEFINE_CONTAINER_PROXY_FACTORY(Globals::DataLocation::Element)

KRATOS_DEFINE_CONTAINER_PROXY_FACTORY(Globals::DataLocation::Condition)

KRATOS_DEFINE_CONTAINER_PROXY_FACTORY(Globals::DataLocation::Constraint)

#undef KRATOS_DEFINE_CONTAINER_PROXY_FACTORY


} // namespace Kratos
