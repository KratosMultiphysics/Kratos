#pragma once

#include "Point.h"

#include <type_traits>

namespace ANurbs {
namespace Internals {

template <typename T>
struct Dimension;

template <typename TScalar, int TDimension>
struct Dimension<Point<TScalar, TDimension>>
{
    static int constexpr value = Point<TScalar, TDimension>::Dimension();
};

template <typename T>
struct Scalar;

template <typename TScalar, int TDimension>
struct Scalar<Point<TScalar, TDimension>>
{
    using type = typename Point<TScalar, TDimension>::ScalarType;
};

template <>
struct Scalar<double>
{
    using type = double;
};

template <typename T>
struct Zero;

template <typename TScalar, int TDimension>
struct Zero<Point<TScalar, TDimension>>
{
    static Point<TScalar, TDimension>
    get()
    {
        return Point<TScalar, TDimension>();
    }
};

template <typename T>
struct Nth
{
    static inline typename Scalar<T>::type
    get(
        const T& vector,
        const int index)
    {
        return vector[index];
    }
};

template <typename T>
struct SquaredNorm
{
    static typename Scalar<T>::type
    get(
        const T& v)
    {
        typename Scalar<T>::type sqSum {0};

        for (int i = 0; i < Dimension<T>::value; i++) {
            sqSum += Nth<T>::get(v, i) * Nth<T>::get(v, i);
        }

        return sqSum;
    }
};

template <typename T, typename = void>
struct Norm
{
    static typename Scalar<T>::type
    get(
        const T& vector)
    {
        typename Scalar<T>::type scale {0};

        for (int i = 0; i < Dimension<T>::value; i++) {
            scale += std::abs(Nth<T>::get(vector, i));
        }

        if (scale == 0) {
            return 0;
        }

        typename Scalar<T>::type scaledSum {0};

        for (int i = 0; i < Dimension<T>::value; i++) {
            const typename Scalar<T>::type si = Nth<T>::get(vector, i) / scale;
            scaledSum += si * si;
        }

        return scale * std::sqrt(scaledSum);
    }
};

template <typename T>
struct Norm<T, typename std::enable_if<std::is_arithmetic<T>::value>::type>
{
    static double
    get(
        const double& scalar)
    {
        return std::abs(scalar);
    }
};

template <typename T, typename = void>
struct Cross;

template <typename T>
struct Cross<T, typename std::enable_if<Dimension<T>::value == 2>::type>
{
    static double
    get(
        const T& u,
        const T& v)
    {
        return Nth<T>::get(v, 0) * Nth<T>::get(u, 1) - Nth<T>::get(v, 1) * Nth<T>::get(u, 0);
    }
};

template <typename T>
struct Cross<T, typename std::enable_if<Dimension<T>::value == 3>::type>
{
    static T
    get(
        const T& u,
        const T& v)
    {
        T result;

        result[0] = Nth<T>::get(v, 1) * Nth<T>::get(u, 2) - Nth<T>::get(v, 2) * Nth<T>::get(u, 1);
        result[1] = Nth<T>::get(v, 2) * Nth<T>::get(u, 0) - Nth<T>::get(v, 0) * Nth<T>::get(u, 2);
        result[2] = Nth<T>::get(v, 0) * Nth<T>::get(u, 1) - Nth<T>::get(v, 1) * Nth<T>::get(u, 0);

        return result;
    }
};

} // namespace Internals

template <typename T>
using ScalarTypeOf = typename Internals::Scalar<T>::type;

template <typename T>
static constexpr int
DimensionOf()
{
    return Internals::Dimension<T>::value;
}

template <typename T>
static ScalarTypeOf<T>
Nth(
    const T& vector,
    const int index)
{
    return Internals::Nth<T>::get(vector, index);
}

template <typename T>
static T
Zero()
{
    return Internals::Zero<T>::get();
}

template <typename T>
static ScalarTypeOf<T>
Dot(
    const T& u,
    const T& v)
{
    ScalarTypeOf<T> result {0};

    for (int i = 0; i < DimensionOf<T>(); i++) {
        result += Nth(u, i) * Nth(v, i);
    }

    return result;
}

template <typename T>
static ScalarTypeOf<T>
SquaredNorm(
    const T& v)
{
    return Internals::SquaredNorm<T>::get(v);
}

template <typename T>
static typename std::result_of<decltype(&Internals::Norm<T>::get)(const T&)>
    ::type
Norm(
    const T& vector)
{
    return Internals::Norm<T>::get(vector);
}

template <typename T>
static typename std::result_of<decltype(&Internals::Cross<T>::get)(const T&,
    const T&)>::type
Cross(
    const T& u,
    const T& v)
{
    return Internals::Cross<T>::get(u, v);
}

} // namespace ANurbs
