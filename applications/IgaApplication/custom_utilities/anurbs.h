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
#include <ANurbs/Integration>

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

namespace Kratos {

template <int TDimension>
using CurveGeometryBase = ANurbs::CurveGeometryBase<
    Kratos::array_1d<double, TDimension>>;

template <int TDimension>
using CurveGeometry = ANurbs::CurveGeometry<
    Kratos::array_1d<double, TDimension>>;

template <int TDimension>
using SurfaceGeometryBase = ANurbs::SurfaceGeometryBase<
    Kratos::array_1d<double, TDimension>>;

template <int TDimension>
using SurfaceGeometry = ANurbs::SurfaceGeometry<
    Kratos::array_1d<double, TDimension>>;

template <int TDimension>
using CurveBase = ANurbs::CurveBase<Kratos::array_1d<double, TDimension>>;

template <int TDimension>
using Curve = ANurbs::Curve<CurveGeometry<TDimension>>;

template <int TDimension>
using SurfaceBase = ANurbs::SurfaceBase<Kratos::array_1d<double, TDimension>>;

template <int TDimension>
using Surface = ANurbs::Surface<SurfaceGeometry<TDimension>>;

template <int TDimension>
using CurveOnSurface = ANurbs::CurveOnSurface<
    Kratos::array_1d<double, 2>, Kratos::array_1d<double, TDimension>>;

template <int TDimension>
using PointOnCurveProjection = ANurbs::PointOnCurveProjection<
    Kratos::array_1d<double, TDimension>>;

using TrimmedSurfaceClipping = ANurbs::TrimmedSurfaceClipping<
    Kratos::array_1d<double, 2>>;

} // namespace Kratos

#endif // !defined(KRATOS_ANURBS_H_INCLUDED)
