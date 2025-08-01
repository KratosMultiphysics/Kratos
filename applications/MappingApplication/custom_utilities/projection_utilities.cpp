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
#include "geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions.h"
#include "geometries/nurbs_surface_geometry.h"

namespace Kratos
{
namespace ProjectionUtilities
{

typedef std::size_t SizeType;
typedef std::size_t IndexType;

typedef Geometry<Node> GeometryType;
typedef typename GeometryType::PointsArrayType PointsArrayType;
typedef PointerVector<Node> ContainerNodeType;
typedef NurbsSurfaceGeometry<3, ContainerNodeType> NurbsSurfaceGeometryType;

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

void FillEquationIdVectorIBRA(const GeometryType::Pointer pGeometry,
                          std::vector<int>& rEquationIds, CoordinatesArrayType rCoordinates)
{
    KRATOS_TRY;
    // Get the polynomial degree of the nurbs surface
    IndexType polynomial_degree_u = pGeometry->PolynomialDegree(0);
    IndexType polynomial_degree_v = pGeometry->PolynomialDegree(1);

    // Get the knot vectors of the nurbs surface and extend them to be consistent 
    std::vector<double> knot_vector_u, knot_vector_v;
    pGeometry->SpansLocalSpace(knot_vector_u, 0);
    pGeometry->SpansLocalSpace(knot_vector_v, 1);
    knot_vector_u.insert(knot_vector_u.begin(), knot_vector_u.front());
    knot_vector_u.insert(knot_vector_u.end(), knot_vector_u.back());
    knot_vector_v.insert(knot_vector_v.begin(), knot_vector_v.front());
    knot_vector_v.insert(knot_vector_v.end(), knot_vector_v.back());

    // shape function container.
    NurbsSurfaceShapeFunction shape_function_container(
        polynomial_degree_u, polynomial_degree_v, 0);
    
    // Transform the knot vectors to the required format for the shape function container
    Vector vector_knot_vector_u(knot_vector_u.size()), vector_knot_vector_v(knot_vector_v.size());   
    for (IndexType i = 0; i < knot_vector_u.size(); ++i) {
        vector_knot_vector_u[i] = knot_vector_u[i];
    }

    for (IndexType i = 0; i < knot_vector_v.size(); ++i) {
        vector_knot_vector_v[i] = knot_vector_v[i];
    }

    shape_function_container.ComputeBSplineShapeFunctionValues(vector_knot_vector_u, vector_knot_vector_v, rCoordinates[0], rCoordinates[1]);

    IndexType num_nonzero_cps = shape_function_container.NumberOfNonzeroControlPoints();

    /// Get List of Control Points
    PointsArrayType nonzero_control_points(num_nonzero_cps);
    std::vector<int> cp_indices = shape_function_container.ControlPointIndices(
        pGeometry->PointsNumberInDirection(0), pGeometry->PointsNumberInDirection(1));
    
    for (IndexType j = 0; j < num_nonzero_cps; j++) {
        KRATOS_DEBUG_ERROR_IF_NOT(pGeometry->pGetPoint(cp_indices[j])->Has(INTERFACE_EQUATION_ID)) << pGeometry->pGetPoint(cp_indices[j]) << " does not have an \"INTERFACE_EQUATION_ID\"" << std::endl;
        rEquationIds.push_back(pGeometry->pGetPoint(cp_indices[j])->GetValue(INTERFACE_EQUATION_ID));
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
        // projection is outside the line, searching the closest point
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

PairingIndex ProjectToIBRA(const GeometryType& rGeometry,
                               const Point& rPointToProject,
                               const double LocalCoordTol,
                               Vector& rShapeFunctionValues,
                               std::vector<int>& rEquationIds,
                               double& rProjectionDistance,
                               const bool ComputeApproximation)
{
    KRATOS_TRY

    // Get the parent geometry of the quadrature point and the geometry type
    const GeometryType& geom_parent = rGeometry.GetGeometryParent(0);
    const auto geom_type = geom_parent.GetGeometryType();

    // Declare and initialize the variables needed for the projection
    CoordinatesArrayType local_curve_coords = ZeroVector(3);
    CoordinatesArrayType local_surface_coords = ZeroVector(3);
    CoordinatesArrayType projected_point_global_coords = ZeroVector(3);
    PairingIndex pairing_index = PairingIndex::Unspecified;

    if (geom_type == GeometryData::KratosGeometryType::Kratos_Brep_Curve_On_Surface || geom_type == GeometryData::KratosGeometryType::Kratos_Nurbs_Curve_On_Surface){
        // Get the nurbs surface geometry
        const GeometryType::Pointer p_nurbs_surface = geom_parent.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

        // Initial value for the non-linear projections step
        std::vector<double> curve_span;
        geom_parent.SpansLocalSpace(curve_span, 0);
        local_curve_coords[0] = (curve_span.front() + curve_span.back()) * 0.5;

        // Try to project the point onto the curve
        if (geom_parent.ProjectionPointGlobalToLocalSpace(rPointToProject, local_curve_coords, 1e-6)){
            pairing_index = PairingIndex::Line_Inside;

            // Provide a proper initial guess for the local coordinates to be given as seed for the non-linear projection step 
            std::vector<double> surface_local_span_u, surface_local_span_v;
            p_nurbs_surface->SpansLocalSpace(surface_local_span_u, 0);
            p_nurbs_surface->SpansLocalSpace(surface_local_span_v, 1);
            local_surface_coords[0] = (surface_local_span_u.front() + surface_local_span_u.back()) * 0.5;
            local_surface_coords[1] = (surface_local_span_v.front() + surface_local_span_v.back()) * 0.5;      

            // Get the global coordinates of the projected point
            geom_parent.GlobalCoordinates(projected_point_global_coords, local_curve_coords);

            // Get the local coordinates of the projected point in the parameter space of the surface
            p_nurbs_surface->ProjectionPointGlobalToLocalSpace(projected_point_global_coords, local_surface_coords, 1e-6);

            // Evaluate the shape functions at the local coordinates and get the equations id vector
            p_nurbs_surface->ShapeFunctionsValues(rShapeFunctionValues, local_surface_coords);
            FillEquationIdVectorIBRA(p_nurbs_surface, rEquationIds, local_surface_coords);
            
            // Get the distance between the projected point and the point to project
            rProjectionDistance = norm_2(rPointToProject - projected_point_global_coords);
        }  else if (!ComputeApproximation) { // If the projection fails and no approximation is allowed, return unspecified
            return pairing_index;
        }
    } else if (geom_type == GeometryData::KratosGeometryType::Kratos_Brep_Surface || geom_type == GeometryData::KratosGeometryType::Kratos_Nurbs_Surface){
        // Get the nurbs surface geometry
        const GeometryType::Pointer p_nurbs_surface = geom_parent.pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

        // Initial value for the non-linear projections step
        std::vector<double> surface_knot_vector_u, surface_knot_vector_v;
        p_nurbs_surface->SpansLocalSpace(surface_knot_vector_u, 0);
        p_nurbs_surface->SpansLocalSpace(surface_knot_vector_v, 1);
        local_surface_coords[0] = (surface_knot_vector_u.front() + surface_knot_vector_u.back()) * 0.5;
        local_surface_coords[1] = (surface_knot_vector_v.front() + surface_knot_vector_v.back()) * 0.5;

        // Try to project the point onto the surface
        if (geom_parent.ProjectionPointGlobalToLocalSpace(rPointToProject, local_surface_coords, 1e-6)){
            pairing_index = PairingIndex::Surface_Inside;

            // Evaluate the shape functions at the local coordinates and get the equations id vector
            p_nurbs_surface->ShapeFunctionsValues(rShapeFunctionValues, local_surface_coords);
            FillEquationIdVectorIBRA(p_nurbs_surface, rEquationIds, local_surface_coords);

            // Get the distance between the projected point and the point to project
            CoordinatesArrayType projected_point_global_coords = ZeroVector(3);
            p_nurbs_surface->GlobalCoordinates(projected_point_global_coords, local_surface_coords);
            rProjectionDistance = norm_2(rPointToProject - projected_point_global_coords);
        } else if (!ComputeApproximation) { // If the projection fails and no approximation is allowed, return unspecified
            return PairingIndex::Unspecified;
        } else if (geom_parent.ProjectionPointGlobalToLocalSpace(rPointToProject, local_surface_coords, LocalCoordTol)) { // If the initial projection fails and an approximation is allowed, try 

            // Evaluate the shape functions at the local coordinates and get the equations id vector
            p_nurbs_surface->ShapeFunctionsValues(rShapeFunctionValues, local_surface_coords);
            FillEquationIdVectorIBRA(p_nurbs_surface, rEquationIds, local_surface_coords);

            // Get the distance between the projected point and the point to project
            CoordinatesArrayType projected_point_global_coords = ZeroVector(3);
            p_nurbs_surface->GlobalCoordinates(projected_point_global_coords, local_surface_coords);
            rProjectionDistance = norm_2(rPointToProject - projected_point_global_coords);
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
    } else if (geom_family == GeometryData::KratosGeometryFamily::Kratos_Quadrature_Geometry) {
        rPairingIndex = ProjectToIBRA(rGeometry, rPointToProject, LocalCoordTol, rShapeFunctionValues, rEquationIds, rProjectionDistance, ComputeApproximation);
        is_full_projection = (rPairingIndex == PairingIndex::Line_Inside  || rPairingIndex == PairingIndex::Surface_Inside);
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
