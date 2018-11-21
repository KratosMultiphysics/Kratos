#pragma once

#include "Polygon.h"
#include "Pointer.h"
#include "CurveTessellation.h"
#include "CurveBase.h"
#include "TrimmedSurfaceClipping.h"

#include <Mapbox/Earcut.hpp>

#include <vector>

namespace mapbox {
namespace util {

template <typename TVector>
struct nth<0, TVector> {
    static inline ANurbs::ScalarTypeOf<TVector>
    get(
        const TVector& vector)
    {
        return ANurbs::Nth(vector, 0);
    };
};

template <typename TVector>
struct nth<1, TVector> {
    static inline ANurbs::ScalarTypeOf<TVector>
    get(
        const TVector& vector)
    {
        return ANurbs::Nth(vector, 1);
    };
};

} // namespace util
} // namespace mapbox

namespace ANurbs {

struct TriangleIndices
{
    int a;
    int b;
    int c;
};

struct QuadIndices
{
    int a;
    int b;
    int c;
    int d;
};

template <typename TVector2>
class PolygonTessellation
{
public:
    using Vector2Type = TVector2;
    using ScalarType = ScalarTypeOf<Vector2Type>;

    using PolygonType = Polygon<Vector2Type>;

private:
    std::vector<int> m_triangles;
    std::vector<int> m_quads;

public:
    PolygonTessellation()
    {
    }

    void
    Compute(
        const PolygonType& polygon)
    {
        m_triangles.clear();
        m_quads.clear();

        using Path = std::vector<Vector2Type>;
        using Paths = std::vector<Path>;

        Paths contours;

        contours.push_back(polygon.outer_path);

        for (const auto& inner_path : polygon.inner_paths) {
            contours.push_back(inner_path);
        }

        // FIXME: Check for quad

        m_triangles = mapbox::earcut<int>(contours);
    }

    int
    NbTriangles() const
    {
        return static_cast<int>(m_triangles.size() / 3);
    }

    TriangleIndices
    Triangle(
        int index) const
    {
        TriangleIndices triangle;

        triangle.a = m_triangles.at(index * 3 + 0);
        triangle.b = m_triangles.at(index * 3 + 1);
        triangle.c = m_triangles.at(index * 3 + 2);

        return triangle;
    }

    int
    NbQuads() const
    {
        return static_cast<int>(m_quads.size() / 4);
    }

    QuadIndices
    Quad(
        int index) const
    {
        QuadIndices quad;

        quad.a = m_quads.at(index * 4 + 0);
        quad.b = m_quads.at(index * 4 + 1);
        quad.c = m_quads.at(index * 4 + 2);
        quad.d = m_quads.at(index * 4 + 3);

        return quad;
    }
};

} // namespace ANurbs