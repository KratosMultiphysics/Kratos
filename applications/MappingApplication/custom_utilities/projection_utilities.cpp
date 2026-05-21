//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "projection_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos
{
namespace ProjectionUtilities
{

typedef std::size_t SizeType;
typedef std::size_t IndexType;

typedef Geometry<Node> GeometryType;

namespace {

void FillEquationIdVector(const GeometryType& rGeometry,
                          std::vector<int>& rEquationIds)
{
    KRATOS_TRY

    const SizeType num_points = rGeometry.PointsNumber();

    if (rEquationIds.size() != num_points) {
        rEquationIds.resize(num_points);
    }

    IndexType point_index = 0;
    for (const auto& r_node : rGeometry.Points()) {
        KRATOS_DEBUG_ERROR_IF_NOT(r_node.Has(INTERFACE_EQUATION_ID)) << r_node << " does not have an \"INTERFACE_EQUATION_ID\"" << std::endl;
        rEquationIds[point_index++] = r_node.GetValue(INTERFACE_EQUATION_ID);
    }

    KRATOS_CATCH("")
}

bool IsBetterProjection(const PairingIndex CurrentPairingIndex,
                        const PairingIndex PairingIndexToCheck,
                        const double CurrentMinDistance,
                        const double DistanceToCheck)
{
    // the projection is better if either:
    // - the new pairing-index is better (i.e. larger) than the current pairing-index
    // - the new pairing-index the same as the current pairing-index but the projection distance is smaller
    return (PairingIndexToCheck > CurrentPairingIndex || (PairingIndexToCheck == CurrentPairingIndex && DistanceToCheck < CurrentMinDistance));
}

}

PairingIndex ProjectOnLine(const GeometryType& rGeometry,
                           const Point& rPointToProject,
                           const double LocalCoordTol,
                           Vector& rShapeFunctionValues,
                           std::vector<int>& rEquationIds,
                           double& rProjectionDistance,
                           const bool ComputeApproximation)
{
    KRATOS_TRY

    Point projected_point;

    rProjectionDistance = std::abs(GeometricalProjectionUtilities::FastProjectOnLine(rGeometry, rPointToProject, projected_point));

    array_1d<double, 3> local_coords;
    PairingIndex pairing_index;

    if (rGeometry.IsInside(projected_point, local_coords, 1e-14)) {
        pairing_index = PairingIndex::Line_Inside;
        rGeometry.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
        FillEquationIdVector(rGeometry, rEquationIds);

    } else if (!ComputeApproximation) {
        return PairingIndex::Unspecified;

    } else if (rGeometry.IsInside(projected_point, local_coords, LocalCoordTol)) {
        pairing_index = PairingIndex::Line_Outside;
        rGeometry.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
        FillEquationIdVector(rGeometry, rEquationIds);

    } else {
        // projection is ouside the line, searching the closest point
        pairing_index = PairingIndex::Closest_Point;
        const double dist_1 = rPointToProject.Distance(rGeometry[0]);
        const double dist_2 = rPointToProject.Distance(rGeometry[1]);

        rEquationIds.resize(1);
        if (dist_1 < dist_2) {
            KRATOS_DEBUG_ERROR_IF_NOT(rGeometry[0].Has(INTERFACE_EQUATION_ID)) << rGeometry[0] << " does not have an \"INTERFACE_EQUATION_ID\"" << std::endl;
            rEquationIds[0] = rGeometry[0].GetValue(INTERFACE_EQUATION_ID);
            rProjectionDistance = dist_1;
        } else {
            KRATOS_DEBUG_ERROR_IF_NOT(rGeometry[1].Has(INTERFACE_EQUATION_ID)) << rGeometry[1] << " does not have an \"INTERFACE_EQUATION_ID\"" << std::endl;
            rEquationIds[0] = rGeometry[1].GetValue(INTERFACE_EQUATION_ID);
            rProjectionDistance = dist_2;
        }

        rShapeFunctionValues.resize(1);
        rShapeFunctionValues[0] = 1.0;
    }

    return pairing_index;

    KRATOS_CATCH("")
}

PairingIndex ProjectOnSurface(const GeometryType& rGeometry,
                     const Point& rPointToProject,
                     const double LocalCoordTol,
                     Vector& rShapeFunctionValues,
                     std::vector<int>& rEquationIds,
                     double& rProjectionDistance,
                     const bool ComputeApproximation)
{
    KRATOS_TRY

    Point projected_point;

    rProjectionDistance = std::abs(GeometricalProjectionUtilities::FastProjectOnGeometry(rGeometry, rPointToProject, projected_point));

    array_1d<double, 3> local_coords;
    PairingIndex pairing_index;

    if (rGeometry.IsInside(projected_point, local_coords, 1e-14)) {
        pairing_index = PairingIndex::Surface_Inside;
        rGeometry.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
        FillEquationIdVector(rGeometry, rEquationIds);

    } else if (!ComputeApproximation) {
        return PairingIndex::Unspecified;

    } else if (rGeometry.IsInside(projected_point, local_coords, LocalCoordTol)) {
        pairing_index = PairingIndex::Surface_Outside;
        rGeometry.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
        FillEquationIdVector(rGeometry, rEquationIds);

    } else { // inter-/extrapolation failed, trying to project on "subgeometries"
        pairing_index = PairingIndex::Unspecified;
        std::vector<int> edge_eq_ids;
        Vector edge_sf_values;
        double edge_distance;

        for (const auto& r_edge : rGeometry.GenerateEdges()) {

            const PairingIndex edge_index = ProjectOnLine(r_edge, rPointToProject, LocalCoordTol, edge_sf_values, edge_eq_ids, edge_distance);

            // check if the current edge gives a better result
            if (IsBetterProjection(pairing_index, edge_index, rProjectionDistance, edge_distance)) {
                pairing_index = edge_index;
                rShapeFunctionValues = edge_sf_values;
                rProjectionDistance = edge_distance;
                rEquationIds = edge_eq_ids;
            }
        }
    }

    return pairing_index;

    KRATOS_CATCH("")
}

PairingIndex ProjectIntoVolume(const GeometryType& rGeometry,
                               const Point& rPointToProject,
                               const double LocalCoordTol,
                               Vector& rShapeFunctionValues,
                               std::vector<int>& rEquationIds,
                               double& rProjectionDistance,
                               const bool ComputeApproximation)
{
    KRATOS_TRY

    array_1d<double, 3> local_coords;
    PairingIndex pairing_index;

    if (rGeometry.IsInside(rPointToProject, local_coords, 1e-14)) {
        pairing_index = PairingIndex::Volume_Inside;
        rGeometry.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
        FillEquationIdVector(rGeometry, rEquationIds);

        rProjectionDistance = rPointToProject.Distance(rGeometry.Center());
        rProjectionDistance /= rGeometry.Volume(); // Normalize Distance by Volume

    } else if (!ComputeApproximation) {
        return PairingIndex::Unspecified;

    } else if (rGeometry.IsInside(rPointToProject, local_coords, LocalCoordTol)) {
        pairing_index = PairingIndex::Volume_Outside;
        rGeometry.ShapeFunctionsValues(rShapeFunctionValues, local_coords);
        FillEquationIdVector(rGeometry, rEquationIds);

        rProjectionDistance = rPointToProject.Distance(rGeometry.Center());
        rProjectionDistance /= rGeometry.Volume(); // Normalize Distance by Volume

    } else { // inter-/extrapolation failed, trying to project on "subgeometries"
        pairing_index = PairingIndex::Unspecified;
        std::vector<int> face_eq_ids;
        Vector face_sf_values;
        double face_distance;

        for (const auto& r_face : rGeometry.GenerateFaces()) {

            const PairingIndex face_index = ProjectOnSurface(r_face, rPointToProject, LocalCoordTol, face_sf_values, face_eq_ids, face_distance);

            // check if the current edge gives a better result
            if (IsBetterProjection(pairing_index, face_index, rProjectionDistance, face_distance)) {
                pairing_index = face_index;
                rShapeFunctionValues = face_sf_values;
                rProjectionDistance = face_distance;
                rEquationIds = face_eq_ids;
            }
        }
    }

    return pairing_index;

    KRATOS_CATCH("")
}

bool ComputeProjection(const GeometryType& rGeometry,
                       const Point& rPointToProject,
                       const double LocalCoordTol,
                       Vector& rShapeFunctionValues,
                       std::vector<int>& rEquationIds,
                       double& rProjectionDistance,
                       PairingIndex& rPairingIndex,
                       const bool ComputeApproximation)
{
    KRATOS_TRY

    const SizeType num_points = rGeometry.PointsNumber();
    const auto geom_family = rGeometry.GetGeometryFamily();
    bool is_full_projection = false;

    if (geom_family == GeometryData::KratosGeometryFamily::Kratos_Linear && num_points == 2) { // linear line
        rPairingIndex = ProjectOnLine(rGeometry, rPointToProject, LocalCoordTol, rShapeFunctionValues, rEquationIds, rProjectionDistance, ComputeApproximation);
        is_full_projection = (rPairingIndex == PairingIndex::Line_Inside);

    } else if ((geom_family == GeometryData::KratosGeometryFamily::Kratos_Triangle      && num_points == 3) || // linear triangle
               (geom_family == GeometryData::KratosGeometryFamily::Kratos_Quadrilateral && num_points == 4)) { // linear quad
        rPairingIndex = ProjectOnSurface(rGeometry, rPointToProject, LocalCoordTol, rShapeFunctionValues, rEquationIds, rProjectionDistance, ComputeApproximation);
        is_full_projection = (rPairingIndex == PairingIndex::Surface_Inside);

    } else if (geom_family == GeometryData::KratosGeometryFamily::Kratos_Tetrahedra ||
               geom_family == GeometryData::KratosGeometryFamily::Kratos_Prism ||
               geom_family == GeometryData::KratosGeometryFamily::Kratos_Pyramid ||
               geom_family == GeometryData::KratosGeometryFamily::Kratos_Hexahedra) { // Volume projection
        rPairingIndex = ProjectIntoVolume(rGeometry, rPointToProject, LocalCoordTol, rShapeFunctionValues, rEquationIds, rProjectionDistance, ComputeApproximation);
        is_full_projection = (rPairingIndex == PairingIndex::Volume_Inside);

    } else if (ComputeApproximation) {
        KRATOS_WARNING_ONCE("Mapper") << "Unsupported type of geometry for projection, trying to use an approximation (Nearest Neighbor)" << std::endl;

        if (rShapeFunctionValues.size() != 1) {
            rShapeFunctionValues.resize(1);
        }
        rShapeFunctionValues[0] = 1.0;

        if (rEquationIds.size() != 1) {
            rEquationIds.resize(1);
        }

        rProjectionDistance = std::numeric_limits<double>::max();
        rPairingIndex = PairingIndex::Closest_Point;
        for (const auto& r_node : rGeometry.Points()) {
            const double dist = rPointToProject.Distance(r_node);
            if (dist < rProjectionDistance) {
                rProjectionDistance = dist;
                KRATOS_DEBUG_ERROR_IF_NOT(r_node.Has(INTERFACE_EQUATION_ID)) << r_node << " does not have an \"INTERFACE_EQUATION_ID\"" << std::endl;
                rEquationIds[0] = r_node.GetValue(INTERFACE_EQUATION_ID);
            }
        }
    }

    KRATOS_DEBUG_ERROR_IF(rPairingIndex != PairingIndex::Unspecified && rShapeFunctionValues.size() != rEquationIds.size()) << "Number of equation-ids is not the same as the number of ShapeFunction values, something went wrong!" << std::endl;

    return is_full_projection;

    KRATOS_CATCH("")
}

} // namespace ProjectionUtilities
} // namespace Kratos.
