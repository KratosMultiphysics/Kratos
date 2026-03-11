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

#include "element_utilities.hpp"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/line_2d_2.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "node_utilities.h"

namespace Kratos
{

Matrix LinearNodalExtrapolator::CalculateElementExtrapolationMatrix(const GeometryType& rGeometry,
                                                                    const Geo::IntegrationPointVectorType& rIntegrationPoints) const
{
    CheckIfGeometryIsSupported(rGeometry);

    const auto  p_lower_order_geometry = CreateLowerOrderGeometry(rGeometry);
    const auto& r_corner_geometry = p_lower_order_geometry ? *p_lower_order_geometry : rGeometry;

    auto result = CalculateExtrapolationMatrixForCornerNodes(rGeometry, rIntegrationPoints, r_corner_geometry);

    if (p_lower_order_geometry) AddRowsForMidsideNodes(rGeometry, result);

    return result;
}

Matrix LinearNodalExtrapolator::CalculateExtrapolationMatrixForCornerNodes(const GeometryType& rGeometry,
                                                                           const Geo::IntegrationPointVectorType& rIntegrationPoints,
                                                                           const GeometryType& rCornerGeometry)
{
    const auto determinants_of_jacobian =
        GeoElementUtilities::EvaluateDeterminantsOfJacobiansAtIntegrationPoints(rIntegrationPoints, rGeometry);

    const auto shape_functions_values_at_integration_points =
        GeoElementUtilities::EvaluateShapeFunctionsAtIntegrationPoints(rIntegrationPoints, rCornerGeometry);

    const auto number_of_corner_nodes = rCornerGeometry.PointsNumber();
    Matrix     quasi_mass_mat         = ZeroMatrix(number_of_corner_nodes, number_of_corner_nodes);
    const auto number_of_integration_points = rIntegrationPoints.size();
    Matrix     node_coefficients(number_of_corner_nodes, number_of_integration_points);
    for (auto i = std::size_t{0}; i < number_of_integration_points; ++i) {
        const auto& N = shape_functions_values_at_integration_points[i];
        quasi_mass_mat += outer_prod(N, N) * determinants_of_jacobian[i] * rIntegrationPoints[i].Weight();
        column(node_coefficients, i) = N * determinants_of_jacobian[i] * rIntegrationPoints[i].Weight();
    }

    double metric_determinant;
    Matrix quasi_mass_mat_inverse;
    MathUtils<>::InvertMatrix(quasi_mass_mat, quasi_mass_mat_inverse, metric_determinant, -1.);

    return prod(quasi_mass_mat_inverse, node_coefficients);
}

void LinearNodalExtrapolator::AddRowsForMidsideNodes(const GeometryType& rGeometry, Matrix& rExtrapolationMatrix)
{
    const auto global_to_local_mapping = NodeUtilities::CreateGlobalToLocalNodeIndexMap(rGeometry.Points());

    rExtrapolationMatrix.resize(rGeometry.PointsNumber(), rExtrapolationMatrix.size2());
    for (const auto& edge : rGeometry.GenerateEdges()) {
        row(rExtrapolationMatrix, global_to_local_mapping.at(edge[2].GetId())) =
            0.5 * (row(rExtrapolationMatrix, global_to_local_mapping.at(edge[0].GetId())) +
                   row(rExtrapolationMatrix, global_to_local_mapping.at(edge[1].GetId())));
    }
}

std::unique_ptr<LinearNodalExtrapolator::GeometryType> LinearNodalExtrapolator::CreateLowerOrderGeometry(const GeometryType& rGeometry)
{
    // Creating lower order geometries is only supported for quadratic geometries.
    if (rGeometry.GetGeometryOrderType() != GeometryData::Kratos_Quadratic_Order) return nullptr;

    switch (rGeometry.GetGeometryFamily()) {
        using enum GeometryData::KratosGeometryFamily;
    case Kratos_Linear:
        return std::make_unique<Line2D2<Node>>(rGeometry(0), rGeometry(1));
    case Kratos_Triangle:
        return std::make_unique<Triangle2D3<Node>>(rGeometry(0), rGeometry(1), rGeometry(2));
    case Kratos_Quadrilateral:
        return std::make_unique<Quadrilateral2D4<Node>>(rGeometry(0), rGeometry(1), rGeometry(2),
                                                        rGeometry(3));
    case Kratos_Tetrahedra:
        return std::make_unique<Tetrahedra3D4<Node>>(rGeometry(0), rGeometry(1), rGeometry(2), rGeometry(3));
    case Kratos_Hexahedra:
        return std::make_unique<Hexahedra3D8<Node>>(rGeometry(0), rGeometry(1), rGeometry(2),
                                                    rGeometry(3), rGeometry(4), rGeometry(5),
                                                    rGeometry(6), rGeometry(7));
    default:
        KRATOS_ERROR << "Cannot create lower order geometry: unsupported family type\n";
    }
}

void LinearNodalExtrapolator::CheckIfGeometryIsSupported(const GeometryType& rGeometry)
{
    KRATOS_ERROR_IF(rGeometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Linear &&
                    rGeometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Triangle &&
                    rGeometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Quadrilateral &&
                    rGeometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Tetrahedra &&
                    rGeometry.GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Hexahedra);

    KRATOS_ERROR_IF(rGeometry.GetGeometryOrderType() != GeometryData::KratosGeometryOrderType::Kratos_Linear_Order &&
                    rGeometry.GetGeometryOrderType() != GeometryData::KratosGeometryOrderType::Kratos_Quadratic_Order);
}

} // namespace Kratos
