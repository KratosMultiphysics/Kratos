#pragma once

#include "SurfaceBase.h"
#include "Pointer.h"
#include "VectorTraits.h"

namespace ANurbs {

template <typename TVector>
class PointOnSurfaceProjection
{
public:
    using VectorType = TVector;
    using ScalarType = ScalarTypeOf<VectorType>;

    using SurfaceBaseType = SurfaceBase<VectorType>;

private:
    struct ParameterPoint
    {
        ScalarType parameterU;
        ScalarType parameterV;
        VectorType point;
    };

private:
    Pointer<SurfaceBaseType> m_surface;
    ParameterPoint m_closestPoint;
    std::vector<ParameterPoint> m_tessellation;
    ScalarType m_tolerance;
    ScalarType m_distance;
    int m_gridU;
    int m_gridV;

public:
    PointOnSurfaceProjection(
        Pointer<SurfaceBaseType> surface)
    : m_surface(surface)
    {
        std::vector<ScalarType> us;

        for (const auto& span : surface->SpansU()) {
            if (span.Length() < 1e-7) {
                continue;
            }

            const int n = surface->DegreeU() + 1;

            for (int i = 0; i < n; i++) {
                const ScalarType u = span.ParameterAtNormalized(1.0 / n);
                us.push_back(u);
            }
        }

        us.push_back(surface->DomainU().T1());

        std::vector<ScalarType> vs;

        for (const auto& span : surface->SpansV()) {
            if (span.Length() < 1e-7) {
                continue;
            }

            const int n = surface->DegreeV() + 1;

            for (int i = 0; i < n; i++) {
                const ScalarType v = span.ParameterAtNormalized(1.0 / n);
                vs.push_back(v);
            }
        }

        vs.push_back(surface->DomainV().T1());

        m_tessellation.reserve(us.size() * vs.size());

        for (const auto u : us) {
            for (const auto v : vs) {
                const auto point = m_surface->PointAt(u, v);
                m_tessellation.push_back({u, v, point});
            }
        }

        m_gridU = static_cast<int>(us.size()) - 1;
        m_gridV = static_cast<int>(vs.size()) - 1;
    }

    Pointer<SurfaceBaseType>
    Surface() const
    {
        return m_surface;
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
    ParameterU() const
    {
        return m_closestPoint.parameterU;
    }

    ScalarType
    ParameterV() const
    {
        return m_closestPoint.parameterV;
    }

    VectorType
    Point() const
    {
        return m_closestPoint.point;
    }

    ScalarType
    Distance() const
    {
        return m_distance;
    }

    void
    Compute(
        const VectorType& sample)
    {
        int minIndex = -1;
        ScalarType minDistance = std::numeric_limits<ScalarType>::max();

        for (int i = 0; i < m_tessellation.size(); i++) {
            const auto pt = m_tessellation[i].point;

            const VectorType v = sample - pt;
            const auto distance = SquaredNorm(v);

            if (distance < minDistance) {
                minDistance = distance;
                minIndex = i;
            }
        }

        // ---

        const int p = minIndex % (m_gridV + 1);
        const int o = minIndex / (m_gridV + 1);

        m_closestPoint = m_tessellation[minIndex];

        if (o != m_gridU && p != m_gridV) {
            const auto pt = TriangleProjection(sample, minIndex, minIndex + 1, minIndex + m_gridV + 1);

            const VectorType v = sample - pt.point;
            const auto distance = SquaredNorm(v);

            if (distance < minDistance) {
                m_closestPoint = pt;
                minDistance = distance;
            }
        }

        if (o != m_gridU && p != 0) {
            const auto pt = TriangleProjection(sample, minIndex, minIndex - 1, minIndex + m_gridV + 1);

            const VectorType v = sample - pt.point;
            const auto distance = SquaredNorm(v);

            if (distance < minDistance) {
                m_closestPoint = pt;
                minDistance = distance;
            }
        }

        if (o != 0 && p != m_gridV) {
            const auto pt = TriangleProjection(sample, minIndex, minIndex + 1, minIndex - m_gridV - 1);

            const VectorType v = sample - pt.point;
            const auto distance = SquaredNorm(v);

            if (distance < minDistance) {
                m_closestPoint = pt;
                minDistance = distance;
            }
        }

        if (o != 0 && p != 0) {
            const auto pt = TriangleProjection(sample, minIndex, minIndex - 1, minIndex - m_gridV - 1);

            const VectorType v = sample - pt.point;
            const auto distance = SquaredNorm(v);

            if (distance < minDistance) {
                m_closestPoint = pt;
                minDistance = distance;
            }
        }

        // ---

        m_closestPoint = Newton(sample, m_closestPoint.parameterU, m_closestPoint.parameterV);
    }

    ParameterPoint
    Newton(
        const VectorType point,
        const ScalarType u,
        const ScalarType v
    )
    {
        ScalarType cu = u;
        ScalarType cv = v;

        int maxits = 5;
        ScalarType eps1 = 0.00001;
        ScalarType eps2 = 0.000005;
        
        ScalarType minu = m_surface->DomainU().T0();
        ScalarType maxu = m_surface->DomainU().T1();
        ScalarType minv = m_surface->DomainV().T0();
        ScalarType maxv = m_surface->DomainV().T1();

        for (int i = 0; i < maxits; i++) {
            const auto s = m_surface->DerivativesAt(cu, cv, 2);

            const VectorType dif = s[0] - point;

            const ScalarType c1v = Norm(dif);
            
            if (c1v < eps1) {
                return {cu, cv, s[0]};
            }

            ScalarType s1_l = Norm(s[1]);
            ScalarType s2_l = Norm(s[2]);
            
            ScalarType c2an = std::abs(Dot(s[1], dif));
            ScalarType c2ad = s1_l * c1v;

            ScalarType c2bn = std::abs(Dot(s[2], dif));
            ScalarType c2bd = s2_l * c1v;

            ScalarType c2av = c2an / c2ad;
            ScalarType c2bv = c2bn / c2bd;

            if (c2av < eps2 && c2bv < eps2) {
                return {cu, cv, s[0]};
            }
            
            ScalarType f = Dot(s[1], dif);
            ScalarType g = Dot(s[2], dif);

            ScalarType J00 = Dot(s[1], s[1]) + Dot(s[3], dif);
            ScalarType J01 = Dot(s[1], s[2]) + Dot(s[4], dif);
            ScalarType J11 = Dot(s[2], s[2]) + Dot(s[5], dif);

            ScalarType det = J00 * J11 - J01 * J01;
            
            ScalarType du = (g * J01 - f * J11) / det;
            ScalarType dv = (f * J01 - g * J00) / det;
            
            cu += du;
            cv += dv;
        }

        return {cu, cv, m_surface->PointAt(cu, cv)};
    }

    ParameterPoint
    TriangleProjection(
        const VectorType point,
        const int& indexA,
        const int& indexB,
        const int& indexC)
    {
        const auto a = m_tessellation[indexA];
        const auto b = m_tessellation[indexB];
        const auto c = m_tessellation[indexC];

        const VectorType u = b.point - a.point;
        const VectorType v = c.point - a.point;
        const VectorType n = Cross(u, v);
        const VectorType w = point - a.point;

        const ScalarType gam = Dot(Cross(u, w), n) / SquaredNorm(n);
        const ScalarType bet = Dot(Cross(w, v), n) / SquaredNorm(n);
        const ScalarType alp = 1.0 - gam - bet;

        ParameterPoint cp;

        cp.parameterU = alp * a.parameterU + bet * b.parameterU +
            gam * c.parameterU;
        cp.parameterV = alp * a.parameterV + bet * b.parameterV +
            gam * c.parameterV;

        cp.point = m_surface->PointAt(cp.parameterU, cp.parameterV);

        return cp;
    }
};

using PointOnSurfaceProjection2D = PointOnSurfaceProjection<Point2D>;
using PointOnSurfaceProjection3D = PointOnSurfaceProjection<Point3D>;

} // namespace ANurbs
