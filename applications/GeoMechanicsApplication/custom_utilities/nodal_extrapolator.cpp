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

    // Sofar this works for 3, 4, 6 and 8 node planar elements
    // for 2 and 3 node line elements the extension is straightforward.
    // for volume elements ( hexa, tetra, wedge ) the midside node interpolation step is more elaborate
    auto p_low_order_geometry = CreateLowerOrderGeometry(rGeometry);
    GeometryType* low_order_geometry = p_low_order_geometry ? p_low_order_geometry.get() : &rGeometry;

    // calculate extrapolation matrix towards corner nodes
    SizeType number_of_low_order_nodes = low_order_geometry->PointsNumber();

    Vector     N              = ZeroVector(number_of_low_order_nodes);
    Matrix     quasi_mass_mat = ZeroMatrix(number_of_low_order_nodes, number_of_low_order_nodes);
    const auto integration_points           = rGeometry.IntegrationPoints(IntegrationMethod);
    const auto number_of_integration_points = integration_points.size();
    Matrix     node_coefficient(number_of_low_order_nodes, number_of_integration_points);
    Vector     determinants_of_jacobian;
    determinants_of_jacobian = rGeometry.DeterminantOfJacobian(determinants_of_jacobian, IntegrationMethod);
    for (IndexType i_gauss_point = 0; i_gauss_point < number_of_integration_points; ++i_gauss_point) {
        // local_coordinates --> isoparametric coordinates or for triangles area coordinates
        const array_1d<double, 3>& r_local_coordinates = integration_points[i_gauss_point].Coordinates();
        // shape function for this i.p.
        low_order_geometry->ShapeFunctionsValues(N, r_local_coordinates);
        quasi_mass_mat += outer_prod(N, N) * determinants_of_jacobian[i_gauss_point] *
                          integration_points[i_gauss_point].Weight();
        column(node_coefficient, i_gauss_point) =
            N * determinants_of_jacobian[i_gauss_point] * integration_points[i_gauss_point].Weight();
    }

    double MetricDet;
    Matrix quasi_mass_mat_inverse;
    MathUtils<double>::InvertMatrix(quasi_mass_mat, quasi_mass_mat_inverse, MetricDet, -1.);

    auto extrapolation_matrix = Matrix{prod(quasi_mass_mat_inverse, node_coefficient)};

    // add rows for midside nodes if needed ( mean value of corner nodes )
    if (rGeometry.PointsNumber() > low_order_geometry->PointsNumber()) {
        extrapolation_matrix.resize(rGeometry.PointsNumber(), extrapolation_matrix.size2());
        std::size_t n_filled = low_order_geometry->PointsNumber();
        for (std::size_t i_row = 0; i_row < rGeometry.PointsNumber() - n_filled; ++i_row) {
            row(extrapolation_matrix, n_filled + i_row) =
                0.5 * (row(extrapolation_matrix, i_row) + row(extrapolation_matrix, (i_row + 1) % n_filled));
        }
    }

    return extrapolation_matrix;
}

std::unique_ptr<NodalExtrapolator::GeometryType> NodalExtrapolator::CreateLowerOrderGeometry(GeometryType& rGeometry) const
{
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
