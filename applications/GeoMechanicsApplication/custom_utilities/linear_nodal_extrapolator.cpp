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
//                   Aron Noordam
//

#include "linear_nodal_extrapolator.h"

#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"

namespace Kratos
{

Matrix LinearNodalExtrapolator::CalculateElementExtrapolationMatrix(const GeometryType& rGeometry,
                                                                    const GeometryData::IntegrationMethod& rIntegrationMethod) const
{
    CheckIfGeometryIsSupported(rGeometry);

    const auto p_lower_order_geometry = CreateLowerOrderGeometry(rGeometry);
    const GeometryType& p_corner_geometry = p_lower_order_geometry ? *p_lower_order_geometry : rGeometry;

    Matrix extrapolation_matrix =
        CalculateExtrapolationMatrixForCornerNodes(rGeometry, rIntegrationMethod, p_corner_geometry);

    if (p_lower_order_geometry) {
        AddRowsForMidsideNodes(rGeometry, extrapolation_matrix);
    }

    return extrapolation_matrix;
}

Matrix LinearNodalExtrapolator::CalculateExtrapolationMatrixForCornerNodes(const GeometryType& rGeometry,
                                                                           const GeometryData::IntegrationMethod& rIntegrationMethod,
                                                                           const GeometryType& rCornerGeometry)
{
    const SizeType number_of_corner_nodes = rCornerGeometry.PointsNumber();

    Matrix      quasi_mass_mat     = ZeroMatrix(number_of_corner_nodes, number_of_corner_nodes);
    const auto& integration_points = rGeometry.IntegrationPoints(rIntegrationMethod);
    const auto  number_of_integration_points = integration_points.size();
    Matrix      node_coefficients(number_of_corner_nodes, number_of_integration_points);
    Vector      determinants_of_jacobian;
    rGeometry.DeterminantOfJacobian(determinants_of_jacobian, rIntegrationMethod);

    const Matrix& shape_functions_values_at_integration_points =
        rCornerGeometry.ShapeFunctionsValues(rIntegrationMethod);

    for (IndexType i = 0; i < number_of_integration_points; ++i) {
        const Vector N = row(shape_functions_values_at_integration_points, i);
        quasi_mass_mat += outer_prod(N, N) * determinants_of_jacobian[i] * integration_points[i].Weight();
        column(node_coefficients, i) = N * determinants_of_jacobian[i] * integration_points[i].Weight();
    }

    double metric_determinant;
    Matrix quasi_mass_mat_inverse;
    MathUtils<double>::InvertMatrix(quasi_mass_mat, quasi_mass_mat_inverse, metric_determinant, -1.);

    return prod(quasi_mass_mat_inverse, node_coefficients);
}

void LinearNodalExtrapolator::AddRowsForMidsideNodes(const GeometryType& rGeometry, Matrix& rExtrapolationMatrix)
{
    const size_t n_filled = rExtrapolationMatrix.size1();
    rExtrapolationMatrix.resize(rGeometry.PointsNumber(), rExtrapolationMatrix.size2());

    if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Triangle ||
        rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Quadrilateral) {
        for (size_t i_row = 0; i_row < rGeometry.PointsNumber() - n_filled; ++i_row) {
            row(rExtrapolationMatrix, n_filled + i_row) =
                0.5 * (row(rExtrapolationMatrix, i_row) + row(rExtrapolationMatrix, (i_row + 1) % n_filled));
        }
    } else if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Tetrahedra &&
               rGeometry.size() == 10) {
        row(rExtrapolationMatrix, 4) = 0.5 * (row(rExtrapolationMatrix, 0) + row(rExtrapolationMatrix, 1));
        row(rExtrapolationMatrix, 5) = 0.5 * (row(rExtrapolationMatrix, 1) + row(rExtrapolationMatrix, 2));
        row(rExtrapolationMatrix, 6) = 0.5 * (row(rExtrapolationMatrix, 2) + row(rExtrapolationMatrix, 0));
        row(rExtrapolationMatrix, 7) = 0.5 * (row(rExtrapolationMatrix, 0) + row(rExtrapolationMatrix, 3));
        row(rExtrapolationMatrix, 8) = 0.5 * (row(rExtrapolationMatrix, 1) + row(rExtrapolationMatrix, 3));
        row(rExtrapolationMatrix, 9) = 0.5 * (row(rExtrapolationMatrix, 2) + row(rExtrapolationMatrix, 3));
    } else if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Hexahedra &&
               rGeometry.size() == 20) {
        row(rExtrapolationMatrix, 8) = 0.5 * (row(rExtrapolationMatrix, 0) + row(rExtrapolationMatrix, 1));
        row(rExtrapolationMatrix, 9) = 0.5 * (row(rExtrapolationMatrix, 1) + row(rExtrapolationMatrix, 2));
        row(rExtrapolationMatrix, 10) = 0.5 * (row(rExtrapolationMatrix, 2) + row(rExtrapolationMatrix, 3));
        row(rExtrapolationMatrix, 11) = 0.5 * (row(rExtrapolationMatrix, 3) + row(rExtrapolationMatrix, 0));
        row(rExtrapolationMatrix, 12) = 0.5 * (row(rExtrapolationMatrix, 0) + row(rExtrapolationMatrix, 4));
        row(rExtrapolationMatrix, 13) = 0.5 * (row(rExtrapolationMatrix, 1) + row(rExtrapolationMatrix, 5));
        row(rExtrapolationMatrix, 14) = 0.5 * (row(rExtrapolationMatrix, 2) + row(rExtrapolationMatrix, 6));
        row(rExtrapolationMatrix, 15) = 0.5 * (row(rExtrapolationMatrix, 3) + row(rExtrapolationMatrix, 7));
        row(rExtrapolationMatrix, 16) = 0.5 * (row(rExtrapolationMatrix, 4) + row(rExtrapolationMatrix, 5));
        row(rExtrapolationMatrix, 17) = 0.5 * (row(rExtrapolationMatrix, 5) + row(rExtrapolationMatrix, 6));
        row(rExtrapolationMatrix, 18) = 0.5 * (row(rExtrapolationMatrix, 6) + row(rExtrapolationMatrix, 7));
        row(rExtrapolationMatrix, 19) = 0.5 * (row(rExtrapolationMatrix, 7) + row(rExtrapolationMatrix, 4));
    } else {
        KRATOS_ERROR << "Unexpected geometry type for AddRowsForMidsideNodes" << std::endl;
    }
}

std::unique_ptr<LinearNodalExtrapolator::GeometryType> LinearNodalExtrapolator::CreateLowerOrderGeometry(const GeometryType& rGeometry)
{
    // Sofar this works for 3, 4, 6 and 8 node planar elements and 4, 8, 10 and 20 node volume elements
    // for 2 and 3 node line elements the extension is straightforward.
    switch (rGeometry.size()) {
    case 6:
        return std::make_unique<Triangle2D3<Node>>(rGeometry(0), rGeometry(1), rGeometry(2));
    case 8:
        // HexaHedra3D8 also has 8 nodes, this should not create a Quadrilateral2D4
        if (rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Quadrilateral) {
            return std::make_unique<Quadrilateral2D4<Node>>(rGeometry(0), rGeometry(1), rGeometry(2),
                rGeometry(3));
        }
        else {
            return nullptr;
        }
    case 10:
        return std::make_unique<Tetrahedra3D4<Node>>(rGeometry(0), rGeometry(1), rGeometry(2), rGeometry(3));

    case 20:
        return std::make_unique<Hexahedra3D8<Node>>(rGeometry(0), rGeometry(1), rGeometry(2),
                                                    rGeometry(3), rGeometry(4), rGeometry(5),
                                                    rGeometry(6), rGeometry(7));
    default:
        return nullptr;
    }
}

void LinearNodalExtrapolator::CheckIfGeometryIsSupported(const GeometryType& rGeometry)
{
    const auto number_of_nodes = rGeometry.size();
    KRATOS_ERROR_IF(rGeometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Triangle &&
                    rGeometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Quadrilateral &&
                    rGeometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Tetrahedra &&
                    rGeometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Hexahedra);

    KRATOS_ERROR_IF(rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Triangle &&
                    (number_of_nodes != 3 && number_of_nodes != 6));

    KRATOS_ERROR_IF(rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Quadrilateral &&
                    (number_of_nodes != 4 && number_of_nodes != 8));

    KRATOS_ERROR_IF(rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Tetrahedra &&
                    (number_of_nodes != 4 && number_of_nodes != 10));

    KRATOS_ERROR_IF(rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Hexahedra &&
                    (number_of_nodes != 8 && number_of_nodes != 20));
}

} // namespace Kratos
