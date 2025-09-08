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
#include "math_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

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

} // namespace Kratos
