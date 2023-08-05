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
#include "includes/model_part.h" // Mesh::NodesContainerType, Mesh::ElementsContainerType, Mesh::ConditionsContainerType

// System includes
#include <type_traits> // remove_reference_t, is_const_v, is_same_v, decay_t
#include <optional> // optional


namespace Kratos {


template <class TEntityProxy>
class ContainerProxy;



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
            return mpEntity.value()->GetSolutionStepValue(rVariable);
        } else {
            return mpEntity.value()->GetValue(rVariable);
        }
    }

    /// @brief Fetch the value corresponding to the input variable in the wrapped entity.
    template <class TVariable, std::enable_if_t</*this is required for SFINAE*/!std::is_same_v<TVariable,void> && TMutable,bool> = true>
    typename TVariable::Type& GetValue(const TVariable& rVariable)
    {
        if constexpr (TLocation == Globals::DataLocation::NodeHistorical) {
            return mpEntity.value()->GetSolutionStepValue(rVariable);
        } else {
            return mpEntity.value()->GetValue(rVariable);
        }
    }

    /// @brief Overwrite the value corresponding to the input variable in the wrapped entity.
    template <class TVariable, std::enable_if_t</*this is required for SFINAE*/!std::is_same_v<TVariable,void> && TMutable,bool> = true>
    void SetValue(const TVariable& rVariable, const typename TVariable::Type& rValue)
    {
        if constexpr (TLocation == Globals::DataLocation::NodeHistorical) {
            mpEntity.value()->GetSolutionStepValue(rVariable) = rValue;
        } else {
            mpEntity.value()->SetValue(rVariable, rValue);
        }
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
                void // <== invalid fallback type; will throw a compile-time error
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

        Iterator operator--(int) noexcept {Iterator copy(mWrapped); --mWrapped; return *this;}

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



#define KRATOS_DEFINE_ENTITY_PROXY_FACTORY(TEntity)                                         \
    template <Globals::DataLocation TLocation>                                              \
    auto MakeProxy(const TEntity& rEntity) {return EntityProxy<TLocation,false>(rEntity);}  \
    template <Globals::DataLocation TLocation>                                              \
    auto MakeProxy(TEntity& rEntity) {return EntityProxy<TLocation,true>(rEntity);}

KRATOS_DEFINE_ENTITY_PROXY_FACTORY(Node)

KRATOS_DEFINE_ENTITY_PROXY_FACTORY(Element)

KRATOS_DEFINE_ENTITY_PROXY_FACTORY(Condition)

#undef KRATOS_DEFINE_ENTITY_PROXY_FACTORY



template <Globals::DataLocation TLocation>
auto MakeProxy(ModelPart& rModelPart)
{
    using TEntityProxy = EntityProxy<TLocation,true>;
    if constexpr (TLocation == Globals::DataLocation::NodeHistorical || TLocation == Globals::DataLocation::NodeNonHistorical) {
        return ContainerProxy<TEntityProxy>(rModelPart.Nodes().begin(), rModelPart.Nodes().end());
    } else if constexpr (TLocation == Globals::DataLocation::Element) {
        return ContainerProxy<TEntityProxy>(rModelPart.Elements().begin(), rModelPart.Elements().end());
    } else if constexpr (TLocation == Globals::DataLocation::Condition) {
        return ContainerProxy<TEntityProxy>(rModelPart.Conditions().begin(), rModelPart.Conditions().end());
    }
}


template <Globals::DataLocation TLocation>
auto MakeProxy(const ModelPart& rModelPart)
{
    using TEntityProxy = EntityProxy<TLocation,false>;
    if constexpr (TLocation == Globals::DataLocation::NodeHistorical || TLocation == Globals::DataLocation::NodeNonHistorical) {
        return ContainerProxy<TEntityProxy>(rModelPart.Nodes().begin(), rModelPart.Nodes().end());
    } else if constexpr (TLocation == Globals::DataLocation::Element) {
        return ContainerProxy<TEntityProxy>(rModelPart.Elements().begin(), rModelPart.Elements().end());
    } else if constexpr (TLocation == Globals::DataLocation::Condition) {
        return ContainerProxy<TEntityProxy>(rModelPart.Conditions().begin(), rModelPart.Conditions().end());
    }
}


} // namespace Kratos
