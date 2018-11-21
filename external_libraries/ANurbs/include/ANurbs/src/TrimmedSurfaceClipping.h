#pragma once

#include "Polygon.h"
#include "Interval.h"
#include "Grid.h"
#include "Pointer.h"
#include "CurveTessellation.h"
#include "CurveBase.h"
#include "SurfaceBase.h"
#include "SurfaceGeometryBase.h"

#include <ClipperLib/clipper.hpp>

#include <vector>

namespace ANurbs {

enum TrimTypes
{
    Empty,
    Full,
    Trimmed
};

template <typename TVector2>
class TrimmedSurfaceClipping
{
public:
    using Vector2Type = TVector2;
    using ScalarType = ScalarTypeOf<Vector2Type>;
    using IntervalType = Interval<ScalarType>;

    using CurveBaseType = CurveBase<Vector2Type>;

    using PolygonType = Polygon<Vector2Type>;

private:
    ClipperLib::Paths m_paths;
    Grid<TrimTypes> m_spanTrimType;
    Grid<std::vector<PolygonType>> m_spanPolygons;
    ScalarType m_scale;
    ScalarType m_tolerance;
    CurveTessellation<Vector2Type> tessellation;
    std::vector<IntervalType> m_spansU;
    std::vector<IntervalType> m_spansV;

    inline ClipperLib::IntPoint
    ToIntPoint(
        const Vector2Type& point) const
    {
        const int x = (int)(point[0] / m_scale);
        const int y = (int)(point[1] / m_scale);

        return {x, y};
    }

    inline Vector2Type
    FromIntPoint(
        const ClipperLib::IntPoint& point) const
    {
        const ScalarType x = point.X * m_scale;
        const ScalarType y = point.Y * m_scale;

        return {x, y};
    }

private:
    static inline bool
    IsRect(
        ClipperLib::Path rect,
        ClipperLib::Path contour,
        std::size_t a,
        std::size_t b,
        std::size_t c,
        std::size_t d)
    {
        if (contour[a] == rect[0] &&
            contour[b] == rect[1] &&
            contour[c] == rect[2] &&
            contour[d] == rect[3]) {
            return true;
        } else {
            return false;
        }
    }

    static inline bool
    IsRect(
        ClipperLib::Path a,
        ClipperLib::Path b)
    {
        return IsRect(a, b, 0, 1, 2, 3) ||
               IsRect(a, b, 1, 2, 3, 0) ||
               IsRect(a, b, 2, 3, 0, 1) ||
               IsRect(a, b, 3, 0, 1, 2) ||
               IsRect(a, b, 3, 2, 1, 0) ||
               IsRect(a, b, 0, 3, 2, 1) ||
               IsRect(a, b, 1, 0, 3, 2) ||
               IsRect(a, b, 2, 1, 0, 3);
    }

    void
    ComputeSpan(
        const int indexU,
        const int indexV,
        const IntervalType spanU,
        const IntervalType spanV)
    {
        ClipperLib::Paths clip(1);
        ClipperLib::PolyTree polytree;

        clip[0] << ToIntPoint({spanU.T0(), spanV.T0()})
                << ToIntPoint({spanU.T1(), spanV.T0()})
                << ToIntPoint({spanU.T1(), spanV.T1()})
                << ToIntPoint({spanU.T0(), spanV.T1()});

        ClipperLib::Clipper clipper;
        clipper.AddPaths(m_paths, ClipperLib::ptSubject, true);
        clipper.AddPaths(clip, ClipperLib::ptClip, true);
        clipper.Execute(ClipperLib::ctIntersection, polytree,
            ClipperLib::pftNonZero, ClipperLib::pftNonZero);

        ClipperLib::PolyNode* polynode = polytree.GetFirst();

        std::vector<PolygonType> regions;

        if (polytree.Total() == 0) {
            m_spanTrimType.SetValue(indexU, indexV, TrimTypes::Empty);
            m_spanPolygons.SetValue(indexU, indexV, {});
            return;
        }

        if (polytree.Total() == 1) {
            const auto& contour = polynode->Contour;

            if (contour.size() == 4 && IsRect(clip[0], contour)) {
                m_spanTrimType.SetValue(indexU, indexV, TrimTypes::Full);
                m_spanPolygons.SetValue(indexU, indexV, {});
                return;
            }
        }

        while (polynode) {
            if (!polynode->IsHole()) {
                const auto& outer_contour = polynode->Contour;
                const auto& inner_contours = polynode->Childs;

                PolygonType region;

                region.outer_path.reserve(outer_contour.size());

                for (const auto& pt : outer_contour) {
                    region.outer_path.push_back(FromIntPoint(pt));
                }

                for (std::size_t i = 0; i < inner_contours.size(); i++) {
                    typename PolygonType::Path inner_path;

                    inner_path.reserve(inner_contours[i]->Contour.size());

                    for (const auto& pt : inner_contours[i]->Contour) {
                        inner_path.push_back(FromIntPoint(pt));
                    }

                    region.inner_paths.push_back(inner_path);
                }

                regions.push_back(region);
            }

            polynode = polynode->GetNext();
        }

        m_spanTrimType.SetValue(indexU, indexV, TrimTypes::Trimmed);
        m_spanPolygons.SetValue(indexU, indexV, regions);
    }

public:
    TrimmedSurfaceClipping(
        const ScalarType tolerance,
        const ScalarType unit)
    : m_tolerance(tolerance)
    , m_scale(unit)
    , m_spanPolygons(1, 1)
    , m_spanTrimType(1, 1)
    {
    }

    void
    Clear()
    {
        m_paths.clear();
    }

    void
    BeginLoop()
    {
        ClipperLib::Path path;
        m_paths.push_back(path);
    }

    void
    EndLoop()
    {
        ClipperLib::CleanPolygon(m_paths.back());
    }

    void
    AddCurve(
        CurveBaseType& curve)
    {
        ClipperLib::Path& path = m_paths.back();

        tessellation.Compute(curve, m_tolerance);

        for (int i = 0; i < tessellation.NbPoints(); i++) {
            auto pt = ToIntPoint(tessellation.Point(i));

            if (i == 0 && path.size() > 0 && path.back() == pt) {
                continue;
            }

            path.push_back(pt);
        }
    }

    template <typename TVector>
    void Compute(
        SurfaceGeometryBase<TVector>& surface)
    {
        Compute(surface.SpansU(), surface.SpansV());
    }

    void Compute(
        const std::vector<IntervalType>& spansU,
        const std::vector<IntervalType>& spansV)
    {
        m_spansU = spansU;
        m_spansV = spansV;

        m_spanTrimType.Resize(NbSpansU(), NbSpansV());
        m_spanPolygons.Resize(NbSpansU(), NbSpansV());

        for (int v = 0; v < NbSpansV(); v++) {
            for (int u = 0; u < NbSpansU(); u++) {
                ComputeSpan(u, v, SpanU(u), SpanV(v));
            }
        }
    }

    int
    NbSpansU() const
    {
        return static_cast<int>(m_spansU.size());
    }

    int
    NbSpansV() const
    {
        return static_cast<int>(m_spansV.size());
    }

    IntervalType
    SpanU(
        const int index) const
    {
        return m_spansU[index];
    }

    IntervalType
    SpanV(
        const int index) const
    {
        return m_spansV[index];
    }

    TrimTypes
    SpanTrimType(
        int indexU,
        int indexV) const
    {
        return m_spanTrimType(indexU, indexV);
    }

    const std::vector<PolygonType>&
    SpanPolygons(
        int indexU,
        int indexV) const
    {
        return m_spanPolygons(indexU, indexV);
    }
};

} // namespace ANurbs