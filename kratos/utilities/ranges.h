//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//
//

#ifndef KRATOS_RANGES_H
#define KRATOS_RANGES_H

// System includes
#include <iterator>
#include <type_traits>


namespace Kratos {


#define KRATOS_DEFINE_RANGE(ITERATOR_TYPE, CONST)                             \
{                                                                            \
public:                                                                      \
    using value_type = typename std::iterator_traits<TIterator>::value_type; \
    using pointer = typename std::iterator_traits<TIterator>::pointer;       \
    using reference = typename std::iterator_traits<TIterator>::reference;   \
    using size_type = std::size_t;                                           \
                                                                             \
    Range() : mBegin(), mEnd() {}                                            \
    Range(ITERATOR_TYPE begin, ITERATOR_TYPE end)                            \
        : mBegin(begin), mEnd(end) {}                                        \
    Range(Range&& rOther) noexcept = default;                                \
    Range(const Range& rOther) noexcept = default;                           \
                                                                             \
    ITERATOR_TYPE begin() CONST noexcept {return mBegin;}                    \
    ITERATOR_TYPE end() CONST noexcept {return mEnd;}                        \
    size_type size() const noexcept {return std::distance(mBegin, mEnd);}    \
    bool empty() const noexcept {return mBegin == mEnd;}                     \
                                                                             \
private:                                                                     \
    ITERATOR_TYPE mBegin;                                                    \
    ITERATOR_TYPE mEnd;                                                      \
}


namespace Detail {
template <class T>
struct IsConstPointer
{};

template <class T>
struct IsConstPointer<T*>
{static constexpr const bool value = false;};

template <class T>
struct IsConstPointer<const T*>
{static constexpr const bool value = true;};

template <class TIterator>
struct IsConstIterator
{
    static constexpr const bool value = IsConstPointer<typename std::iterator_traits<TIterator>::pointer>::value;
};
} // namespace Detail


/**
 *  @brief Class representing a view into a subrange of a container.
 *  @tparam TIterator Iterator type of the target container.
 */
template <class TIterator, class IsConstRange = void>
class Range
KRATOS_DEFINE_RANGE(TIterator, ); // class Range (non-const version)

/**
 *  @brief Class representing a view into a subrange of an immutable container.
 *  @tparam TIterator Iterator type of the target container.
 */
template <class TIterator>
class Range<TIterator, typename std::enable_if<Detail::IsConstIterator<TIterator>::value>::type>
KRATOS_DEFINE_RANGE(TIterator, const); // class Range (const version)

#undef KRATOS_DEFINE_RANGE


/**
 *  @brief Iterator providing access to the keys of an std::map or std::unordered_map.
 *  @tparam TIterator Iterator type of the map.
 */
template <class TIterator>
class MapKeyIterator
{
public:
    using value_type = typename TIterator::value_type::first_type;

    /// const key_type* if const_iterator else key_type*
    using pointer = typename std::conditional<Detail::IsConstIterator<TIterator>::value,const value_type*,value_type*>::type;

    /// const key_type& if const_iterator else key_type&
    using reference = typename std::conditional<Detail::IsConstIterator<TIterator>::value,const value_type&,value_type&>::type;

    using difference_type = typename TIterator::difference_type;

    using iterator_category = std::forward_iterator_tag;

    MapKeyIterator() = default;

    MapKeyIterator(TIterator Wrapped)
        : mWrapped(Wrapped)
    {}

    MapKeyIterator(MapKeyIterator&& rOther) noexcept = default;

    MapKeyIterator(const MapKeyIterator& rOther) noexcept = default;

    reference operator*() {return mWrapped->first;}

    pointer operator->() {return &mWrapped->first;}

    MapKeyIterator& operator++() {++mWrapped; return *this;}

    MapKeyIterator operator++(int) {MapKeyIterator copy(*this); ++(*this); return copy;}

    friend bool operator==(MapKeyIterator Left, MapKeyIterator Right) {return Left.mWrapped == Right.mWrapped;}

    friend bool operator!=(MapKeyIterator Left, MapKeyIterator Right) {return !(Left == Right);}

private:
    TIterator mWrapped;
}; // class MapKeyIterator


/**
 *  @brief Create a view on the keys of an std::map or std::unordered_map.
 *  @note This is a convenience function to avoid having to specify template parameters.
 */
template <class TMap>
Range<MapKeyIterator<typename TMap::const_iterator>> MakeConstMapKeyRange(const TMap& rMap)
{
    return Range<MapKeyIterator<typename TMap::const_iterator>>(rMap.begin(), rMap.end());
}


} // namespace Kratos

#endif // KRATOS_RANGES_H
