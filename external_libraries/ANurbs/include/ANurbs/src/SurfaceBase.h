#pragma once

#include "Interval.h"
#include "VectorTraits.h"

#include <vector>

namespace ANurbs {

template <typename TVector>
class SurfaceBase
{
public:
    using SurfaceBaseType = SurfaceBase<TVector>;

    using VectorType = TVector;
    using ScalarType = ScalarTypeOf<VectorType>;
    using IntervalType = Interval<ScalarType>;

public:
    static constexpr int
    Dimension()
    {
        return DimensionOf<VectorType>();
    }

    virtual int
    DegreeU() const = 0;

    virtual int
    DegreeV() const = 0;

    virtual IntervalType
    DomainU() const = 0;

    virtual IntervalType
    DomainV() const = 0;

    virtual VectorType
    PointAt(
        const ScalarType u,
        const ScalarType v) const = 0;

    virtual std::vector<VectorType>
    DerivativesAt(
        const ScalarType u,
        const ScalarType v,
        const int order) const = 0;

    virtual std::vector<IntervalType>
    SpansU() = 0;

    virtual std::vector<IntervalType>
    SpansV() = 0;
};

using SurfaceBase1D = SurfaceBase<Point1D>;
using SurfaceBase2D = SurfaceBase<Point2D>;
using SurfaceBase3D = SurfaceBase<Point3D>;

} // namespace ANurbs
