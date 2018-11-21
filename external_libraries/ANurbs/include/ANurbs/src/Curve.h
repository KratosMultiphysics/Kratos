#pragma once

#include "CurveBase.h"
#include "CurveGeometry.h"
#include "Pointer.h"

namespace ANurbs {

template <typename TCurveGeometry>
class Curve
    : public CurveBase<typename TCurveGeometry::VectorType>
{
private:
    using CurveGeometryType = TCurveGeometry;

public:
    using CurveBaseType = CurveBase<typename TCurveGeometry::VectorType>;

    using CurveType = Curve<TCurveGeometry>;

    using typename CurveBaseType::IntervalType;
    using typename CurveBaseType::ScalarType;
    using typename CurveBaseType::VectorType;

private:
    Pointer<CurveGeometryType> m_curveGeometry;
    IntervalType m_domain;

public:
    Curve(
        Pointer<CurveGeometryType> geometry)
        : Curve(geometry, geometry->Domain())
    {
    }

    Curve(
        Pointer<CurveGeometryType> geometry,
        const ScalarType t0,
        const ScalarType t1)
        : m_curveGeometry(geometry)
        , m_domain(t0, t1)
    {
    }

    Curve(
        Pointer<CurveGeometryType> geometry,
        const IntervalType& domain)
        : m_curveGeometry(geometry)
        , m_domain(domain)
    {
    }

    Pointer<CurveGeometryType>
    CurveGeometry() const
    {
        return m_curveGeometry;
    }

    int
    Degree() const override
    {
        return m_curveGeometry->Degree();
    }

    IntervalType
    Domain() const override
    {
        return m_domain;
    }

    VectorType
    PointAt(
        const ScalarType t) const override
    {
        return m_curveGeometry->PointAt(t);
    }

    std::vector<VectorType>
    DerivativesAt(
        const ScalarType t,
        const int order) const override
    {
        return m_curveGeometry->DerivativesAt(t, order);
    }

    std::vector<IntervalType>
    Spans() const override
    {
        auto knots = CurveGeometry()->Knots();

        int firstSpan = Knots::UpperSpan(Degree(), knots, Domain().T0());
        int lastSpan = Knots::LowerSpan(Degree(), knots, Domain().T1());

        int nbSpans = lastSpan - firstSpan + 1;

        std::vector<IntervalType> result(nbSpans);

        for (int i = 0; i < nbSpans; i++) {
            ScalarType t0 = CurveGeometry()->Knot(firstSpan + i);
            ScalarType t1 = CurveGeometry()->Knot(firstSpan + i + 1);

            result[i] = IntervalType(t0, t1);
        }

        return result;
    }
};

using Curve1D = Curve<CurveGeometry1D>;
using Curve2D = Curve<CurveGeometry2D>;
using Curve3D = Curve<CurveGeometry3D>;

} // namespace ANurbs
