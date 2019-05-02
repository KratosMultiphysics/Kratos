#pragma once

#include <functional>
#include <type_traits>

namespace ANurbs {
namespace Util {

template <typename TContainer>
struct CurveWeights;

template <typename TWeight>
struct CurveWeights<std::vector<TWeight>>
{
    using Type = TWeight;
    using ContainerType = std::vector<TWeight>;

    static Type
    Get(
        const ContainerType& container,
        const int index)
    {
        return container[index];
    }
};

template <typename TContainer>
struct CurveWeights
{
    using Type = typename std::result_of<TContainer(int)>::type;
    using ContainerType = TContainer;

    static Type
    Get(
        const ContainerType& container,
        const int index)
    {
        return container(index);
    }
};

template <typename TContainer>
struct SurfaceWeights;

template <typename TContainer>
struct SurfaceWeights
{
    using Type =
        typename std::result_of<TContainer(int, int)>::type;
    using ContainerType = TContainer;

    static Type
    Get(
        const ContainerType& container,
        const int indexU,
        const int indexV)
    {
        return container(indexU, indexV);
    }
};

} // Util
} // ANurbs
