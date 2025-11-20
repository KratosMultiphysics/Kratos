// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#include "geometry_utilities.h"
#include "geometries/geometry_data.h"
#include "math_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;

std::size_t GetNumberOfCornerPoints(GeometryData::KratosGeometryFamily GeometryFamily)
{
    switch (GeometryFamily) {
        using enum GeometryData::KratosGeometryFamily;
    case Kratos_Linear:
        return 2;
    case Kratos_Triangle:
        return 3;
    case Kratos_Quadrilateral:
        return 4;
    default:
        KRATOS_ERROR << "The specified geometry family is not supported for getting the number of "
                        "corner points.\n";
    }
}

std::size_t GetNumberOfEdgePoints(GeometryData::KratosGeometryFamily    GeometryFamily,
                                  GeometryData::KratosGeometryOrderType GeometryOrder)
{
    switch (GeometryOrder) {
        using enum GeometryData::KratosGeometryOrderType;
    case Kratos_Linear_Order:
        return 0;
    case Kratos_Quadratic_Order:
        return GetNumberOfCornerPoints(GeometryFamily);
    case Kratos_Cubic_Order:
        return 2 * GetNumberOfCornerPoints(GeometryFamily);
    case Kratos_Quartic_Order:
        return 3 * GetNumberOfCornerPoints(GeometryFamily);
    default:
        KRATOS_ERROR
            << "The specified geometry order type is not supported for getting the number of "
               "edge points.\n";
    }
}

template <std::random_access_iterator InputIt>
void ReverseNodes(InputIt                               Begin,
                  InputIt                               End,
                  GeometryData::KratosGeometryFamily    GeometryFamily,
                  GeometryData::KratosGeometryOrderType GeometryOrderType)
{
    // For line geometries we want to reverse all 'corner points', while for surfaces we don't
    // change the starting node, but only reverse the order of the rest of the corner points.
    auto begin_of_corner_points =
        GeometryFamily == GeometryData::KratosGeometryFamily::Kratos_Linear ? Begin : Begin + 1;

    const auto number_of_corner_points = GetNumberOfCornerPoints(GeometryFamily);
    KRATOS_ERROR_IF(static_cast<int>(number_of_corner_points) > std::distance(Begin, End))
        << "Number of nodes for reversal is too small for the geometry family and order type "
           "specified.\n";
    auto end_of_corner_points = Begin + number_of_corner_points;
    std::reverse(begin_of_corner_points, end_of_corner_points);

    auto end_of_edge_points = End;
    if (GeometryFamily != GeometryData::KratosGeometryFamily::Kratos_Linear) {
        // For non-line geometries, there could be internal points as well, so we only reverse the
        // edge points here. For line geometries, all remaining points will be edge points.
        const auto number_of_edge_points = GetNumberOfEdgePoints(GeometryFamily, GeometryOrderType);
        KRATOS_ERROR_IF(static_cast<int>(number_of_edge_points) > std::distance(end_of_corner_points, End))
            << "Number of nodes for reversal is too small for the geometry family and order type "
               "specified.\n";
        end_of_edge_points = end_of_corner_points + number_of_edge_points;
    }

    std::reverse(end_of_corner_points, end_of_edge_points);

    if (std::distance(end_of_edge_points, End) > 1) {
        std::reverse(end_of_edge_points + 1, End);
    }
}

} // namespace

namespace Kratos
{

Matrix GeometryUtilities::Calculate2DRotationMatrixForLineGeometry(const Geometry<Node>& rGeometry,
                                                                   const array_1d<double, 3>& rLocalCoordinate)
{
    // Since the shape functions depend on one coordinate only
    // for lines, the jacobian only has one column.
    Matrix jacobian;
    rGeometry.Jacobian(jacobian, rLocalCoordinate);
    const auto tangential_vector = GeoMechanicsMathUtilities::Normalized(Vector{column(jacobian, 0)});

    // clang-format off
    Matrix result(2, 2);
    result <<= tangential_vector[0], -tangential_vector[1],
               tangential_vector[1],  tangential_vector[0];
    // clang-format on

    return result;
}

Matrix GeometryUtilities::Calculate3DRotationMatrixForPlaneGeometry(const Geometry<Node>& rGeometry,
                                                                    const array_1d<double, 3>& rLocalCoordinate)
{
    Matrix jacobian;
    rGeometry.Jacobian(jacobian, rLocalCoordinate);
    const auto tangential_vector_1 = GeoMechanicsMathUtilities::Normalized(Vector{column(jacobian, 0)});
    const auto tangential_vector_2 = GeoMechanicsMathUtilities::Normalized(
        Vector{Vector{column(jacobian, 1)} -
               inner_prod(Vector{column(jacobian, 1)}, tangential_vector_1) * tangential_vector_1});
    Vector normal_vector(3);
    MathUtils<>::CrossProduct(normal_vector, tangential_vector_1, tangential_vector_2);
    normal_vector = GeoMechanicsMathUtilities::Normalized(normal_vector);

    // clang-format off
    Matrix result(3, 3);
    result <<= tangential_vector_1[0], tangential_vector_2[0], normal_vector[0],
               tangential_vector_1[1], tangential_vector_2[1], normal_vector[1],
               tangential_vector_1[2], tangential_vector_2[2], normal_vector[2];
    // clang-format on

    return result;
}

std::vector<std::size_t> GeometryUtilities::GetNodeIdsFromGeometry(const Geometry<Node>& rGeometry)
{
    std::vector<std::size_t> result;
    result.reserve(rGeometry.size());
    std::ranges::transform(rGeometry, std::back_inserter(result),
                           [](const auto& rNode) { return rNode.Id(); });
    return result;
}

void GeometryUtilities::ReverseNodes(PointerVector<Node>&                  rNodes,
                                     GeometryData::KratosGeometryFamily    GeometryFamily,
                                     GeometryData::KratosGeometryOrderType GeometryOrderType)
{
    ::ReverseNodes(rNodes.ptr_begin(), rNodes.ptr_end(), GeometryFamily, GeometryOrderType);
}

void GeometryUtilities::ReverseNodes(std::vector<std::size_t>&             rNodeIds,
                                     GeometryData::KratosGeometryFamily    GeometryFamily,
                                     GeometryData::KratosGeometryOrderType GeometryOrderType)
{
    ::ReverseNodes(rNodeIds.begin(), rNodeIds.end(), GeometryFamily, GeometryOrderType);
}

} // namespace Kratos
