#pragma once

#include "CurveBase.h"
#include "CurveGeometry.h"
#include "CurveTessellation.h"
#include "SurfaceGeometry.h"
#include "Pointer.h"

#include <utility>

namespace ANurbs {

template <typename TVector>
class PointOnCurveProjection
{
public:
    using VectorType = TVector;
    using ScalarType = ScalarTypeOf<VectorType>;

    using CurveBaseType = CurveBase<VectorType>;

private:
    struct ParameterPoint
    {
        ScalarType parameter;
        VectorType point;
    };

private:
    CurveTessellation<VectorType> m_tessellation;
    Pointer<CurveBaseType> m_curve;
    ScalarType m_tessellationFlatness;
    ScalarType m_tolerance;
    ScalarType m_parameter;
    VectorType m_point;

public:
    PointOnCurveProjection(
        Pointer<CurveBaseType> curve,
        const ScalarType& tolerance)
    : m_tessellation()
    , m_curve(curve)
    , m_tessellationFlatness(1e-3)
    , m_tolerance(tolerance)
    {
        // create polyline
        m_tessellation.Compute(*Curve(), TessellationFlatness());
    }

    Pointer<CurveBaseType>
    Curve() const
    {
        return m_curve;
    }

    ScalarType
    TessellationFlatness() const
    {
        return m_tessellationFlatness;
    }

    void
    SetTessellationFlatness(
        ScalarType value)
    {
        m_tessellationFlatness = value;
    }

    ScalarType
    Tolerance() const
    {
        return m_tolerance;
    }

    void
    SetTolerance(
        ScalarType value)
    {
        m_tolerance = value;
    }

    ScalarType
    Parameter() const
    {
        return m_parameter;
    }
    
    VectorType
    Point() const
    {
        return m_point;
    }

    void
    Compute(
        const VectorType& sample)
    {
        auto domain = Curve()->Domain();

        // closest point to polyline

        ScalarType closestParameter;
        VectorType closestPoint;

        ScalarType closestSqrDistance = std::numeric_limits<ScalarType>::max();

        for (int i = 1; i < m_tessellation.NbPoints(); i++) {
            ScalarType t0 = m_tessellation.Parameter(i - 1);
            VectorType point0 = m_tessellation.Point(i - 1);

            ScalarType t1 = m_tessellation.Parameter(i);
            VectorType point1 = m_tessellation.Point(i);

            auto projection = LinearProjection(sample, point0, point1, t0, t1);

            ScalarType parameter = projection.parameter;
            VectorType point = projection.point;

            auto sqrDistance = SquaredNorm(VectorType(point - sample));

            if (sqrDistance < closestSqrDistance) {
                closestSqrDistance = sqrDistance;
                closestParameter = parameter;
                closestPoint = point;
            }
        }

        // newton-raphson

        int maxIterations = 5;
        double eps1 = Tolerance();
        double eps2 = Tolerance() * 5;

        for (int i = 0; i < maxIterations; i++) {
            auto f = Curve()->DerivativesAt(closestParameter, 2);

            VectorType dif = f[0] - sample;

            ScalarType c1v = Norm(dif);

            ScalarType c2n = Dot(f[1], dif);
            ScalarType c2d = Norm(f[1]) * c1v;
            ScalarType c2v = c2d != 0 ? c2n / c2d : 0;

            bool c1 = c1v < eps1;
            bool c2 = std::abs(c2v) < eps2;

            if (c1 || c2) { // FIXME: check if 'or' is correct (NURBS Book P.231)
                break;
            }

            ScalarType delta = Dot(f[1], dif) / (Dot(f[2], dif)
                + SquaredNorm(f[1]));

            ScalarType nextParameter = closestParameter - delta;

            // FIXME: out-of-domain check

            // FIXME: 3. condition: (nextParameter - closestParameter) * f[1].Norm();

            closestParameter = domain.Clamp(nextParameter);
        }

        closestPoint = Curve()->PointAt(closestParameter);
        
        closestSqrDistance = SquaredNorm(VectorType(sample - closestPoint));

        VectorType pointAtT0 = Curve()->PointAt(domain.T0());

        if (SquaredNorm(VectorType(sample - pointAtT0)) < closestSqrDistance) {
            m_parameter = domain.T0();
            m_point = pointAtT0;
            return;
        }

        VectorType pointAtT1 = Curve()->PointAt(domain.T1());

        if (SquaredNorm(VectorType(sample - pointAtT1)) < closestSqrDistance) {
            m_parameter = domain.T1();
            m_point = pointAtT1;
            return;
        }

        m_parameter = closestParameter;
        m_point = closestPoint;
    }

private:
    static ParameterPoint
    LinearProjection(
        const VectorType& point,
        const VectorType& a,
        const VectorType& b,
        const ScalarType& t0,
        const ScalarType& t1
    )
    {
        VectorType dif = b - a;
        ScalarType l = SquaredNorm(dif);

        if (l < 1e-14) {
            return {t0, a};
        }

        VectorType o = a;
        VectorType r = dif * (1.0 / l);
        VectorType o2pt = point - o;
        ScalarType do2ptr = Dot(o2pt, r);

        if (do2ptr < 0) {
            return {t0, a};
        }
        
        if (do2ptr > 1) {
            return {t1, b};
        }

        ScalarType t = t0 + (t1 - t0) * do2ptr / l;
        VectorType closestPoint = o + r * do2ptr;

        return {t, closestPoint};
    }
};

using PointOnCurveProjection2D = PointOnCurveProjection<Point2D>;
using PointOnCurveProjection3D = PointOnCurveProjection<Point3D>;

} // namespace ANurbs
