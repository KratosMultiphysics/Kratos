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

#include "nodal_extrapolator.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/triangle_2d_3.h"

namespace Kratos
{

Matrix NodalExtrapolator::CalculateElementExtrapolationMatrix(GeometryType& rGeometry,
                                                              GeometryData::IntegrationMethod IntegrationMethod) const
{
    CheckIfGeometryIsSupported(rGeometry);

    auto          p_lower_order_geometry = CreateLowerOrderGeometry(rGeometry);
    GeometryType& p_corner_geometry = p_lower_order_geometry ? *p_lower_order_geometry : rGeometry;

    Matrix extrapolation_matrix =
        CalculateExtrapolationMatrixForCornerNodes(rGeometry, IntegrationMethod, p_corner_geometry);

    if (p_lower_order_geometry) {
        AddRowsForMidsideNodes(rGeometry, extrapolation_matrix);
    }

    return extrapolation_matrix;
}

Matrix NodalExtrapolator::CalculateExtrapolationMatrixForCornerNodes(
    const NodalExtrapolator::GeometryType& rGeometry,
    GeometryData::IntegrationMethod&       IntegrationMethod,
    const NodalExtrapolator::GeometryType& rCornerGeometry) const
{
    SizeType number_of_corner_nodes = rCornerGeometry.PointsNumber();

    Matrix      quasi_mass_mat     = ZeroMatrix(number_of_corner_nodes, number_of_corner_nodes);
    const auto& integration_points = rGeometry.IntegrationPoints(IntegrationMethod);
    const auto  number_of_integration_points = integration_points.size();
    Matrix      node_coefficient(number_of_corner_nodes, number_of_integration_points);
    Vector      determinants_of_jacobian;
    determinants_of_jacobian = rGeometry.DeterminantOfJacobian(determinants_of_jacobian, IntegrationMethod);

    const Matrix& shape_functions_values_at_integration_points =
        rCornerGeometry.ShapeFunctionsValues(IntegrationMethod);

    for (IndexType i = 0; i < number_of_integration_points; ++i) {
        const Vector N = row(shape_functions_values_at_integration_points, i);
        quasi_mass_mat += outer_prod(N, N) * determinants_of_jacobian[i] * integration_points[i].Weight();
        column(node_coefficient, i) = N * determinants_of_jacobian[i] * integration_points[i].Weight();
    }

    double metric_determinant;
    Matrix quasi_mass_mat_inverse;
    MathUtils<double>::InvertMatrix(quasi_mass_mat, quasi_mass_mat_inverse, metric_determinant, -1.);

    return Matrix{prod(quasi_mass_mat_inverse, node_coefficient)};
}

void NodalExtrapolator::AddRowsForMidsideNodes(const NodalExtrapolator::GeometryType& rGeometry,
                                               Matrix& extrapolation_matrix) const
{
    const size_t n_filled = extrapolation_matrix.size2();
    extrapolation_matrix.resize(rGeometry.PointsNumber(), extrapolation_matrix.size2());
    for (size_t i_row = 0; i_row < rGeometry.PointsNumber() - n_filled; ++i_row) {
        row(extrapolation_matrix, n_filled + i_row) =
            0.5 * (row(extrapolation_matrix, i_row) + row(extrapolation_matrix, (i_row + 1) % n_filled));
    }
}

std::unique_ptr<NodalExtrapolator::GeometryType> NodalExtrapolator::CreateLowerOrderGeometry(GeometryType& rGeometry) const
{
    // Sofar this works for 3, 4, 6 and 8 node planar elements
    // for 2 and 3 node line elements the extension is straightforward.
    // for volume elements ( hexa, tetra, wedge ) the midside node interpolation step is more elaborate
    switch (rGeometry.size()) {
    case 6:
        return std::make_unique<Triangle2D3<Node>>(rGeometry(0), rGeometry(1), rGeometry(2));
    case 8:
        return std::make_unique<Quadrilateral2D4<Node>>(rGeometry(0), rGeometry(1), rGeometry(2),
                                                        rGeometry(3));
    default:
        return nullptr;
    }
}

void NodalExtrapolator::CheckIfGeometryIsSupported(const GeometryType& r_this_geometry) const
{
    const auto number_of_nodes = r_this_geometry.size();
    KRATOS_ERROR_IF(r_this_geometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Triangle &&
                    r_this_geometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Quadrilateral);

    KRATOS_ERROR_IF(r_this_geometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Triangle &&
                    (number_of_nodes != 3 && number_of_nodes != 6));

    KRATOS_ERROR_IF(r_this_geometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Quadrilateral &&
                    (number_of_nodes != 4 && number_of_nodes != 8));
}

} // namespace Kratos
