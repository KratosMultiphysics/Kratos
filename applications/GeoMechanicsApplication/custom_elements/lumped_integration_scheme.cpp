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

#include <geometries/geometry.h>
#include <geometries/quadrilateral_2d_8.h>
#include <geometries/triangle_2d_6.h>
#include <includes/node.h>

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
        Geometry<Node>::PointsArrayType points;
        points.push_back(Node::Pointer(new Node(1, 0.0, 0.0, 0.0)));
        points.push_back(Node::Pointer(new Node(2, 1.0, 0.0, 0.0)));
        points.push_back(Node::Pointer(new Node(3, 0.0, 1.0, 0.0)));
        points.push_back(Node::Pointer(new Node(4, 0.5, 0.0, 0.0)));
        points.push_back(Node::Pointer(new Node(5, 0.5, 0.5, 0.0)));
        points.push_back(Node::Pointer(new Node(6, 0.0, 0.5, 0.0)));
        auto p_geom          = std::make_shared<Triangle2D6<Node>>(points);
        auto lumping_factors = Vector{};
        p_geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::DIAGONAL_SCALING);
        return {{Point(0.0, 0.0), lumping_factors[0]}, {Point(1.0, 0.0), lumping_factors[1]},
                {Point(0.0, 1.0), lumping_factors[2]}, {Point(0.5, 0.0), lumping_factors[3]},
                {Point(1.0, 0.5), lumping_factors[4]}, {Point(0.0, 0.5), lumping_factors[5]}};
    }
    case 8: {
        Geometry<Node>::PointsArrayType points;
        points.push_back(Node::Pointer(new Node(1, -1.0, -1.0, 0.0)));
        points.push_back(Node::Pointer(new Node(2, 1.0, -1.0, 0.0)));
        points.push_back(Node::Pointer(new Node(3, 1.0, 1.0, 0.0)));
        points.push_back(Node::Pointer(new Node(4, -1.0, 1.0, 0.0)));
        points.push_back(Node::Pointer(new Node(5, 0.0, -1.0, 0.0)));
        points.push_back(Node::Pointer(new Node(6, 1.0, 0.0, 0.0)));
        points.push_back(Node::Pointer(new Node(7, 0.0, 1.0, 0.0)));
        points.push_back(Node::Pointer(new Node(8, -1.0, 0.0, 0.0)));
        auto p_geom          = std::make_shared<Quadrilateral2D8<Node>>(points);
        auto lumping_factors = Vector{};
        p_geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::DIAGONAL_SCALING);

        return {{Point(-1.0, -1.0), lumping_factors[0]}, {Point(1.0, -1.0), lumping_factors[1]},
                {Point(1.0, 1.0), lumping_factors[2]},   {Point(-1.0, 1.0), lumping_factors[3]},
                {Point(0.0, -1.0), lumping_factors[4]},  {Point(1.0, 0.0), lumping_factors[5]},
                {Point(0.0, 1.0), lumping_factors[6]},   {Point(-1.0, 0.0), lumping_factors[7]}};
    }
    default:
        KRATOS_ERROR << "Can't construct Lumped integration scheme: no support for "
                     << NumberOfPoints << " point(s)" << std::endl;
    }
}

} // namespace Kratos
