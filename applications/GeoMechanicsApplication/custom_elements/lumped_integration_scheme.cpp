// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#include "lumped_integration_scheme.h"

#include "geometries/geometry.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/triangle_2d_6.h"
#include "includes/node.h"

namespace Kratos
{

LumpedIntegrationScheme::LumpedIntegrationScheme(std::size_t NumberOfPoints)
    : mIntegrationPoints(CreateIntegrationPoints(NumberOfPoints))
{
}

std::size_t LumpedIntegrationScheme::GetNumberOfIntegrationPoints() const
{
    return mIntegrationPoints.size();
}

const Geo::IntegrationPointVectorType& LumpedIntegrationScheme::GetIntegrationPoints() const
{
    return mIntegrationPoints;
}

Geo::IntegrationPointVectorType LumpedIntegrationScheme::CreateIntegrationPoints(std::size_t NumberOfPoints)
{
    // Integration scheme used for planar interface elements in GeoMechanicsApplication
    // locations coincide with node-pairs, numbering follows the node-pairs too, weights as in lumped mass matrix from diagonal scaling of a surface element ( see geometry.h )
    switch (NumberOfPoints) {
    case 3:
        return {{Point(0.0, 0.0), 1.0 / 3.0}, {Point(1.0, 0.0), 1.0 / 3.0}, {Point(0.0, 1.0), 1.0 / 3.0}};
    case 4:
        return {{Point(-1.0, -1.0), 1.0 / 4.0},
                {Point(1.0, -1.0), 1.0 / 4.0},
                {Point(1.0, 1.0), 1.0 / 4.0},
                {Point(-1.0, 1.0), 1.0 / 4.0}};
    case 6: {
        const auto points = std::vector<Point>{Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0),
                                               Point(0.5, 0.0), Point(0.5, 0.5), Point(0.0, 0.5)};
        std::vector<std::size_t> node_ids(points.size());
        std::iota(node_ids.begin(), node_ids.end(), std::size_t{1});
        auto node_ptrs = Geometry<Node>::PointsArrayType{};
        for (auto i = std::size_t{0}; i < node_ids.size(); ++i) {
            node_ptrs.push_back(make_intrusive<Node>(node_ids[i], points[i]));
        }

        auto geom            = Triangle2D6<Node>(node_ptrs);
        auto lumping_factors = Vector{};
        geom.LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::DIAGONAL_SCALING);

        auto make_integration_point = [](const Point& rPoint, double Weight) {
            return Geo::IntegrationPointType{rPoint, Weight};
        };
        Geo::IntegrationPointVectorType result;
        std::transform(points.cbegin(), points.cend(), lumping_factors.cbegin(),
                       std::back_inserter(result), make_integration_point);
        return result;
    }
    case 8: {
        const auto points = std::vector<Point>{Point(-1.0, -1.0), Point(1.0, -1.0), Point(1.0, 1.0),
                                               Point(-1.0, 1.0),  Point(0.0, -1.0), Point(1.0, 0.0),
                                               Point(0.0, 1.0),   Point(-1.0, 0.0)};
        std::vector<std::size_t> node_ids(points.size());
        std::iota(node_ids.begin(), node_ids.end(), std::size_t{1});
        auto node_ptrs = Geometry<Node>::PointsArrayType{};
        for (auto i = std::size_t{0}; i < node_ids.size(); ++i) {
            node_ptrs.push_back(make_intrusive<Node>(node_ids[i], points[i]));
        }

        auto geom            = Quadrilateral2D8<Node>(node_ptrs);
        auto lumping_factors = Vector{};
        geom.LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::DIAGONAL_SCALING);

        auto make_integration_point = [](const Point& rPoint, double Weight) {
            return Geo::IntegrationPointType{rPoint, Weight};
        };
        Geo::IntegrationPointVectorType result;
        std::transform(points.cbegin(), points.cend(), lumping_factors.cbegin(),
                       std::back_inserter(result), make_integration_point);
        return result;
    }
    default:
        KRATOS_ERROR << "Can't construct Lumped integration scheme: no support for "
                     << NumberOfPoints << " point(s)" << std::endl;
    }
}

} // namespace Kratos
