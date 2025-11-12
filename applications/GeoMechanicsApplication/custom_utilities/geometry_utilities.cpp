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
[[nodiscard]] std::size_t GetNumberOfCornerPoints(const GeometryData::KratosGeometryFamily& rGeometryFamily)
{
    switch (rGeometryFamily) {
        using enum GeometryData::KratosGeometryFamily;
    case Kratos_Linear:
        return 2;
    case Kratos_Triangle:
        return 3;
    case Kratos_Quadrilateral:
        return 4;
    default:
        KRATOS_ERROR << "The geometry family of the mid-geometry is not supported\n";
    }
}
}


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

std::vector<std::size_t> GeometryUtilities::GetReversedNodeIdsForGeometryFamily(
    const GeometryData::KratosGeometryFamily& rGeometryFamily, const std::vector<std::size_t>& rInitialNodeIds)
{
    auto result = rInitialNodeIds;
    std::reverse(result.begin(), result.begin() + GetNumberOfCornerPoints(rGeometryFamily));
    return result;
}

} // namespace Kratos
