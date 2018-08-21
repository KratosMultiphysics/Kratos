/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

#if !defined(KRATOS_ANURBS_H_INCLUDED)
#define KRATOS_ANURBS_H_INCLUDED

// System includes

// External includes
#include <ANurbs/Core>

// Project includes

namespace ANurbs {
namespace Internals {

template <typename TScalar, int TDimension>
struct Dimension<Kratos::array_1d<TScalar, TDimension>>
{
    static int constexpr value = TDimension;
};

template <typename TScalar, int TDimension>
struct Scalar<Kratos::array_1d<TScalar, TDimension>>
{
    using type = TScalar;
};

template <typename TScalar, int TDimension>
struct Zero<Kratos::array_1d<TScalar, TDimension>>
{
    static Kratos::array_1d<TScalar, TDimension> get()
    {
        Kratos::array_1d<TScalar, TDimension> zero;
        for (std::size_t i = 0; i < TDimension; i++) {
            zero[i] = 0;
        }
        return zero;
    }
};

} // namespace Internals
} // namespace ANurbs

#endif // !defined(KRATOS_ANURBS_H_INCLUDED)
