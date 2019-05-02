#pragma once

#include "Pointer.h"
#include "SurfaceBase.h"
#include "SurfaceGeometry.h"

namespace ANurbs {

template <typename TSurfaceGeometry>
class Surface
    : public SurfaceBase<typename TSurfaceGeometry::VectorType>
{
private:
    using SurfaceGeometryType = TSurfaceGeometry;

public:
    using SurfaceBaseType = SurfaceBase<typename TSurfaceGeometry::VectorType>;

    using SurfaceType = Surface<TSurfaceGeometry>;

    using typename SurfaceBaseType::IntervalType;
    using typename SurfaceBaseType::ScalarType;
    using typename SurfaceBaseType::VectorType;

private:
    Pointer<SurfaceGeometryType> m_surfaceGeometry;
    IntervalType m_domainU;
    IntervalType m_domainV;

public:
    Surface(
        Pointer<SurfaceGeometryType> geometry)
        : Surface(geometry, geometry->DomainU(), geometry->DomainV())
    {
    }

    Surface(
        Pointer<SurfaceGeometryType> geometry,
        const IntervalType& domainU,
        const IntervalType& domainV)
        : m_surfaceGeometry(geometry)
        , m_domainU(domainU)
        , m_domainV(domainV)
    {
    }

    Pointer<SurfaceGeometryType>
    SurfaceGeometry() const
    {
        return m_surfaceGeometry;
    }

    int
    DegreeU() const override
    {
        return m_surfaceGeometry->DegreeU();
    }

    int
    DegreeV() const override
    {
        return m_surfaceGeometry->DegreeV();
    }

    IntervalType
    DomainU() const override
    {
        return m_domainU;
    }

    IntervalType
    DomainV() const override
    {
        return m_domainV;
    }

    VectorType
    PointAt(
        const ScalarType u,
        const ScalarType v) const override
    {
        return m_surfaceGeometry->PointAt(u, v);
    }

    std::vector<VectorType>
    DerivativesAt(
        const ScalarType u,
        const ScalarType v,
        const int order) const override
    {
        return m_surfaceGeometry->DerivativesAt(u, v, order);
    }

    std::vector<IntervalType>
    SpansU() override
    {
        auto knots = m_surfaceGeometry->KnotsU();

        int firstSpan = Knots::UpperSpan(DegreeU(), knots, DomainU().T0());
        int lastSpan = Knots::LowerSpan(DegreeU(), knots, DomainU().T1());

        int nbSpans = lastSpan - firstSpan + 1;

        std::vector<IntervalType> result(nbSpans);

        for (int i = 0; i < nbSpans; i++) {
            ScalarType u0 = SurfaceGeometry()->KnotU(firstSpan + i);
            ScalarType u1 = SurfaceGeometry()->KnotU(firstSpan + i + 1);

            result[i] = IntervalType(u0, u1);
        }

        return result;
    }

    std::vector<IntervalType>
    SpansV() override
    {
        auto knots = m_surfaceGeometry->KnotsV();

        int firstSpan = Knots::UpperSpan(DegreeV(), knots, DomainV().T0());
        int lastSpan = Knots::LowerSpan(DegreeV(), knots, DomainV().T1());

        int nbSpans = lastSpan - firstSpan + 1;

        std::vector<IntervalType> result(nbSpans);

        for (int i = 0; i < nbSpans; i++) {
            ScalarType v0 = SurfaceGeometry()->KnotV(firstSpan + i);
            ScalarType v1 = SurfaceGeometry()->KnotV(firstSpan + i + 1);

            result[i] = IntervalType(v0, v1);
        }

        return result;
    }
};

using Surface1D = Surface<SurfaceGeometry1D>;
using Surface2D = Surface<SurfaceGeometry2D>;
using Surface3D = Surface<SurfaceGeometry3D>;

} // namespace ANurbs
