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

} // namespace Kratos
