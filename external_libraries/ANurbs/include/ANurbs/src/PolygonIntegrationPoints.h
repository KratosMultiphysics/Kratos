#pragma once

#include "IntegrationPoints.h"
#include "Interval.h"
#include "Polygon.h"
#include "PolygonTessellation.h"

#include <stdexcept>
#include <vector>

namespace ANurbs {

template <typename TVector2>
class PolygonIntegrationPoints
{
public:
    using Vector2Type = TVector2;
    using ScalarType = ScalarTypeOf<Vector2Type>;

    using IntegrationPoint2Type = IntegrationPoint2<ScalarType>;
    using IntegrationPointBarycentricType = IntegrationPointBarycentric<ScalarType>;

private:
    std::vector<IntegrationPoint2Type> m_integrationPoints;

    PolygonTessellation<Vector2Type> m_tessellation;

public:
    void
    Compute(
        const size_t degree,
        const Polygon<Vector2Type>& polygon)
    {
        m_tessellation.Compute(polygon);

        const auto& xiaogimb =
            IntegrationPoints<ScalarType>::XiaoGimbutas(degree);

        const int nbIntegrationPoints =
            m_tessellation.NbTriangles() * xiaogimb.size();

        m_integrationPoints.resize(nbIntegrationPoints);

        auto integrationPoint = m_integrationPoints.begin();

        for (int i = 0; i < m_tessellation.NbTriangles(); i++) {
            const auto triangle = m_tessellation.Triangle(i);

            const Vector2Type vertexA = polygon.Vertex(triangle.a);
            const Vector2Type vertexB = polygon.Vertex(triangle.b);
            const Vector2Type vertexC = polygon.Vertex(triangle.c);

            const Vector2Type vectorAB = vertexB - vertexA;
            const Vector2Type vectorAC = vertexC - vertexA;

            const ScalarType area = 0.5 * Norm(Cross(vectorAB, vectorAC));

            for (const auto& normalizedPoint : xiaogimb) {
                const auto uv = vertexA * normalizedPoint.a +
                    vertexB * normalizedPoint.b + vertexC * normalizedPoint.c;

                integrationPoint->u = uv[0];
                integrationPoint->v = uv[1];
                integrationPoint->weight = area * normalizedPoint.weight;
                integrationPoint++;
            }
        }
    }

    int
    NbIntegrationPoints() const
    {
        return static_cast<int>(m_integrationPoints.size());
    }

    IntegrationPoint2Type
    IntegrationPoint(
        const int index) const
    {
        return m_integrationPoints[index];
    }
};

} // namespace ANurbs
