//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

// System includes
#include <cmath>
#include <limits>

// Project includes
#include "custom_utilities/sbm_boundary_integration_utility.h"
#include "includes/variables.h"
#include "integration/line_gauss_legendre_integration_points.h"
#include "utilities/math_utils.h"

namespace Kratos {

SbmBoundaryIntegrationUtility::IntegrationMode
SbmBoundaryIntegrationUtility::GetIntegrationMode(
    const Condition& rCondition)
{
    KRATOS_ERROR_IF_NOT(rCondition.Has(SURROGATE_BOUNDARY_PROJECTION))
        << rCondition.Info() << " #" << rCondition.Id()
        << " is missing SURROGATE_BOUNDARY_PROJECTION." << std::endl;

    const Matrix& r_projection_data = rCondition.GetValue(
        SURROGATE_BOUNDARY_PROJECTION);
    KRATOS_ERROR_IF(r_projection_data.size1() != NumberOfIntegrationPoints)
        << rCondition.Info() << " #" << rCondition.Id()
        << " expects one SURROGATE_BOUNDARY_PROJECTION row for each of the "
        << NumberOfIntegrationPoints << " boundary integration points."
        << std::endl;

    if (r_projection_data.size2() == ShiftedProjectionColumns) {
        return IntegrationMode::ShiftedSurrogate;
    }
    if (r_projection_data.size2() == ReconstructedQuadratureColumns) {
        return IntegrationMode::ReconstructedBoundary;
    }

    KRATOS_ERROR
        << rCondition.Info() << " #" << rCondition.Id()
        << " expects SURROGATE_BOUNDARY_PROJECTION with "
        << ShiftedProjectionColumns << " columns (classical shifted SBM) or "
        << ReconstructedQuadratureColumns
        << " columns (reconstructed-boundary Gap-SBM)." << std::endl;
    return IntegrationMode::ShiftedSurrogate;
}

void SbmBoundaryIntegrationUtility::CheckConditionData(
    const Condition& rCondition,
    const bool RequireReconstructedBoundary)
{
    KRATOS_ERROR_IF_NOT(rCondition.Has(SURROGATE_BOUNDARY_FACE_COORDINATES))
        << rCondition.Info() << " #" << rCondition.Id()
        << " is missing SURROGATE_BOUNDARY_FACE_COORDINATES." << std::endl;
    const Matrix& r_face_coordinates = rCondition.GetValue(
        SURROGATE_BOUNDARY_FACE_COORDINATES);
    KRATOS_ERROR_IF(
        r_face_coordinates.size1() != 2 ||
        r_face_coordinates.size2() != 3)
        << rCondition.Info() << " #" << rCondition.Id()
        << " expects a 2x3 surrogate-face coordinate matrix." << std::endl;

    const IntegrationMode mode = GetIntegrationMode(rCondition);
    KRATOS_ERROR_IF(
        RequireReconstructedBoundary &&
        mode != IntegrationMode::ReconstructedBoundary)
        << rCondition.Info() << " #" << rCondition.Id()
        << " requires reconstructed-boundary Gap-SBM quadrature." << std::endl;
}

SbmBoundaryIntegrationUtility::IntegrationPointData
SbmBoundaryIntegrationUtility::CalculateIntegrationPointData(
    const Condition& rCondition,
    const IndexType IntegrationPointIndex)
{
    KRATOS_ERROR_IF(IntegrationPointIndex >= NumberOfIntegrationPoints)
        << rCondition.Info() << " #" << rCondition.Id()
        << " requested invalid boundary integration point "
        << IntegrationPointIndex << "." << std::endl;

    const IntegrationMode mode = GetIntegrationMode(rCondition);
    const Matrix& r_face_coordinates = rCondition.GetValue(
        SURROGATE_BOUNDARY_FACE_COORDINATES);
    const Matrix& r_projection_data = rCondition.GetValue(
        SURROGATE_BOUNDARY_PROJECTION);
    const auto& r_integration_points =
        LineGaussLegendreIntegrationPoints2::IntegrationPoints();

    const double tangent_x =
        r_face_coordinates(1, 0) - r_face_coordinates(0, 0);
    const double tangent_y =
        r_face_coordinates(1, 1) - r_face_coordinates(0, 1);
    const double surrogate_face_length = std::hypot(tangent_x, tangent_y);
    KRATOS_ERROR_IF(
        surrogate_face_length <= std::numeric_limits<double>::epsilon())
        << rCondition.Info() << " #" << rCondition.Id()
        << " has a zero-length surrogate face." << std::endl;

    const double local_coordinate =
        r_integration_points[IntegrationPointIndex].X();
    const double first_weight = 0.5 * (1.0 - local_coordinate);
    const double second_weight = 0.5 * (1.0 + local_coordinate);
    Point surrogate_point(0.0, 0.0, 0.0);
    for (IndexType d = 0; d < 3; ++d) {
        surrogate_point[d] = first_weight * r_face_coordinates(0, d) +
            second_weight * r_face_coordinates(1, d);
    }

    IntegrationPointData data;
    data.Mode = mode;
    data.EvaluationPoint = surrogate_point;
    data.Normal[0] = tangent_y / surrogate_face_length;
    data.Normal[1] = -tangent_x / surrogate_face_length;
    data.Weight = 0.5 * surrogate_face_length *
        r_integration_points[IntegrationPointIndex].Weight();

    if (mode == IntegrationMode::ReconstructedBoundary) {
        // Gap-SBM already provides a complete physical-boundary quadrature
        // row: coordinates, oriented unit normal and physical line weight.
        for (IndexType d = 0; d < 3; ++d) {
            data.EvaluationPoint[d] =
                r_projection_data(IntegrationPointIndex, d);
        }
        data.Normal[0] = r_projection_data(IntegrationPointIndex, 3);
        data.Normal[1] = r_projection_data(IntegrationPointIndex, 4);
        data.Weight = r_projection_data(IntegrationPointIndex, 6);
        KRATOS_ERROR_IF(data.Weight < 0.0)
            << rCondition.Info() << " #" << rCondition.Id()
            << " has a negative reconstructed-boundary integration weight."
            << std::endl;
    }

    CalculateParentShapeFunctions(
        rCondition,
        data.EvaluationPoint,
        data.ShapeFunctions,
        data.ShapeFunctionGradients);
    data.TraceShapeFunctions = data.ShapeFunctions;

    if (mode == IntegrationMode::ShiftedSurrogate) {
        // The classical SBM trace is the first-order Taylor expansion from
        // the surrogate Gauss point to its true-boundary projection.
        const array_1d<double, 3> shift{
            r_projection_data(IntegrationPointIndex, 0) - surrogate_point.X(),
            r_projection_data(IntegrationPointIndex, 1) - surrogate_point.Y(),
            r_projection_data(IntegrationPointIndex, 2) - surrogate_point.Z()};
        for (IndexType i = 0; i < data.TraceShapeFunctions.size(); ++i) {
            data.TraceShapeFunctions[i] +=
                data.ShapeFunctionGradients(i, 0) * shift[0] +
                data.ShapeFunctionGradients(i, 1) * shift[1];
        }
    }

    return data;
}

void SbmBoundaryIntegrationUtility::CalculateParentShapeFunctions(
    const Condition& rCondition,
    const Point& rGlobalPoint,
    Vector& rN,
    Matrix& rDN_DX)
{
    const auto& r_geometry = rCondition.GetGeometry();
    GeometryType::CoordinatesArrayType local_coordinates = ZeroVector(3);
    r_geometry.PointLocalCoordinates(local_coordinates, rGlobalPoint);
    r_geometry.ShapeFunctionsValues(rN, local_coordinates);

    Matrix DN_De;
    r_geometry.ShapeFunctionsLocalGradients(DN_De, local_coordinates);
    Matrix jacobian;
    r_geometry.Jacobian(jacobian, local_coordinates);
    Matrix inverse_jacobian;
    double determinant_jacobian = 0.0;
    MathUtils<double>::InvertMatrix(
        jacobian, inverse_jacobian, determinant_jacobian);
    KRATOS_ERROR_IF(
        std::abs(determinant_jacobian) <=
        std::numeric_limits<double>::epsilon())
        << rCondition.Info() << " #" << rCondition.Id()
        << " has a singular owner element." << std::endl;

    rDN_DX.resize(DN_De.size1(), inverse_jacobian.size2(), false);
    noalias(rDN_DX) = prod(DN_De, inverse_jacobian);
}

} // namespace Kratos
