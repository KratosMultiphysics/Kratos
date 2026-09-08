//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti
//

// System includes
#include <fstream>
#include <iomanip>

// External includes

// Project includes
#include "iga_mapping_intersection_utilities.h"
#include "custom_utilities/mapping_triangulation_utilities.h"

// Data structures for doing spatial search 
#include "utilities/function_parser_utility.h"
#include "spatial_containers/bins_dynamic.h"

namespace Kratos
{
namespace
{

using CoordinatesArrayType =
    IgaMappingIntersectionUtilities::CoordinatesArrayType;
using GeometryType = IgaMappingIntersectionUtilities::GeometryType;
using PatchCacheMap = IgaMappingIntersectionUtilities::PatchCacheMap;
using NurbsSurfaceType = NurbsSurfaceGeometry<3, PointerVector<Node>>;

constexpr double parametric_point_tolerance = 1e-8;

bool AddUniquePoint(
    std::vector<CoordinatesArrayType>& rPoints,
    const CoordinatesArrayType& rPoint);

bool IsInsideInitialFEMTriangle(
    const CoordinatesArrayType& rPoint,
    const GeometryType& rFEMGeometry);

void AddEnclosedDomainCorners(
    const NurbsSurfaceType& rNurbsSurface,
    const GeometryType& rFEMGeometry,
    std::vector<CoordinatesArrayType>& rPoints);

bool RecoverUnprojectedTriangleIntersection(
    const GeometryType& rMasterGeometry,
    const GeometryType& rFEMGeometry,
    const NurbsSurfaceType& rNurbsSurface,
    const PatchCacheMap& rPatchCache,
    const double SearchRadius,
    std::vector<CoordinatesArrayType>& rPolygon);

// Writes the three vertices of a parametric triangle and its BRep patch ID.
template<class TTriangleType>
void WriteParametricTriangleToTextFile(
    std::ostream& rOutput,
    const TTriangleType& rTriangle,
    const std::size_t BrepId);

// Maps a parametric triangle to physical space and writes its vertices.
template<class TTriangleType>
void WritePhysicalTriangleToTextFile(
    std::ostream& rOutput,
    const TTriangleType& rTriangle,
    const GeometryType& rBrepGeometry);

} // unnamed namespace

// Pairs each FEM curve condition with its closest projectable IGA curve. The
// resulting pairs are stored as coupling geometries in the result model part.
void IgaMappingIntersectionUtilities::CreateIgaFEMCouplingGeometriesOnCurve(
    ModelPart& rModelPartDomainA,
    ModelPart& rModelPartDomainB,
    const bool& rIsOriginIga,
    ModelPart& rModelPartResult,
    double Tolerance)
{
    const auto& r_iga_model_part = rIsOriginIga ? rModelPartDomainA : rModelPartDomainB;
    const auto& r_fem_model_part = rIsOriginIga ? rModelPartDomainB : rModelPartDomainA;

    for (auto& fem_cond : r_fem_model_part.Conditions()) {
        const auto& line_geom = fem_cond.GetGeometry();
        const Point line_center = line_geom.Center();
        
        double min_distance = std::numeric_limits<double>::max();
        Condition::Pointer p_closest_iga_condition = nullptr;

        // Get the closest condition on the IGA interface
        for (auto& iga_cond : r_iga_model_part.Conditions()) {
            const auto& brep_curve_on_surface_geom = iga_cond.GetGeometry();

            CoordinatesArrayType local_coords = ZeroVector(3);
            CoordinatesArrayType projected_point = ZeroVector(3);

            bool success = brep_curve_on_surface_geom.ProjectionPointGlobalToLocalSpace(line_center, local_coords, Tolerance);
            if (!success) continue;

            brep_curve_on_surface_geom.GlobalCoordinates(projected_point, local_coords);
            double distance = norm_2(line_center - projected_point);

            if (distance < min_distance) {
                min_distance = distance;
                p_closest_iga_condition = r_iga_model_part.pGetCondition(iga_cond.Id());
            }
        }

        // Create a coupling geometry relating both sides
        if (p_closest_iga_condition) {
            rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                p_closest_iga_condition->pGetGeometry(),
                fem_cond.pGetGeometry()
            ));
        }
    }
}

// Creates integration-point coupling conditions on previously paired IGA and FEM curves.
void IgaMappingIntersectionUtilities::CreateIgaFEMQuadraturePointsOnCurve(
    ModelPart& rModelPartCoupling,
    double Tolerance)
{
   const ModelPart &rParentModelPart = rModelPartCoupling.GetParentModelPart();

   // Loop over the coupling geometries and for each one create the integration points on the interface
   for (auto geometry_itr = rModelPartCoupling.GeometriesBegin();
             geometry_itr != rModelPartCoupling.GeometriesEnd();
             ++geometry_itr)
    {
        IntegrationPointsArrayType integration_points;
        IntegrationInfo integration_info = geometry_itr->GetDefaultIntegrationInfo();

        geometry_itr->CreateIntegrationPoints(integration_points, integration_info);

        // Vector of coupling geometries which stores the origin and destination quadrature point geometries
        GeometriesArrayType origin_and_destination_quadrature_points_geometries(integration_points.size());
        IndexType number_of_shape_functions_derivatives = 3;

        if (integration_points.size() != 0){
            geometry_itr->CreateQuadraturePointGeometries(origin_and_destination_quadrature_points_geometries, number_of_shape_functions_derivatives, integration_points, integration_info);

            const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
                                    ? 1
                                    : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;

            for (IndexType i = 0; i < origin_and_destination_quadrature_points_geometries.size(); ++i)
            {
                rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                    id + i, origin_and_destination_quadrature_points_geometries(i)));
            }
        }
    }
}

// Creates candidate coupling geometries between FEM surface conditions and nearby
// IGA patches. A spatial cache limits the patches considered for each FEM condition.
void IgaMappingIntersectionUtilities::CreateIgaFEMCouplingGeometriesOnSurface(
        ModelPart &rModelPartDomainA,
        ModelPart &rModelPartDomainB,
        ModelPart &rModelPartResult,
        bool is_origin_iga,
        const double search_radius,
        PatchCacheMap& rPatchCache)
{
    // Get a pointer to the brep surface starting from the underlying element or condition (quadrature point)
    std::vector<IndexType> patches_id;

    if (is_origin_iga == true){
        std::unordered_set<IndexType> patch_ids;

        // Looping over the elements to get the patch ids 
        for (const auto& r_elem : rModelPartDomainA.Elements())
        {
            const auto& r_geom = r_elem.GetGeometry();
            const auto& r_parent = r_geom.GetGeometryParent(0);
            patch_ids.insert(r_parent.Id());
        }
        
        // Looping over the conditions to get the patch ids
        for (const auto& r_cond : rModelPartDomainA.Conditions())
        {
            const auto& r_geom = r_cond.GetGeometry();
            const auto& r_parent = r_geom.GetGeometryParent(0);
            patch_ids.insert(r_parent.Id());
        }

        patches_id.assign(patch_ids.begin(), patch_ids.end());
    }

    const IndexType patch_divisions = 100;
    rPatchCache = IgaMappingIntersectionUtilities::BuildPatchCaches(
        patches_id, rModelPartDomainA, patch_divisions);

    /* - Iterate over the elements
        - For each element, we should decide which patches are the ones that the element has a projection on 
        - Create a new list called patchs_id_projection depending on the element considered */
    for (auto condition_b_itr = rModelPartDomainB.ConditionsBegin(); condition_b_itr != rModelPartDomainB.ConditionsEnd(); ++condition_b_itr){
        CoordinatesArrayType element_center = condition_b_itr->GetGeometry().Center();

        std::vector<IndexType> patches_with_probable_projection_id = IgaMappingIntersectionUtilities::GetPatchesWithProbableProjection(patches_id, rPatchCache, element_center, search_radius);

        for (IndexType i = 0; i < patches_with_probable_projection_id.size(); i++){
            // Get a pointer to the brep surface geometry
            auto p_brep_surface = rModelPartDomainA.GetRootModelPart().pGetGeometry(patches_with_probable_projection_id[i]);

            // Cast the nurbs surface pointer and brep surface to the derived class
            auto brep_surface_cast = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(p_brep_surface);

            // Creating one coupling geometry for each combination of load condition and patch surface
            rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                    p_brep_surface, condition_b_itr->pGetGeometry()));
        }       
    }
}

// Samples each IGA patch and builds spatial bins for fast projection searches.
// Returns a cache containing the sampled points, bins, and parametric-grid information.
IgaMappingIntersectionUtilities::PatchCacheMap IgaMappingIntersectionUtilities::BuildPatchCaches(
    const std::vector<IndexType>& patches_id,
    const ModelPart& rModelPartIga,
    const IndexType n_div)
{
    PatchCacheMap cache;
    cache.reserve(patches_id.size());

    for (const IndexType patch_id : patches_id) {
        PatchSearchCache pc;

        // Get BrepSurface
        auto p_geom = rModelPartIga.pGetGeometry(patch_id);

        auto p_brep = dynamic_pointer_cast<
            BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>
        >(p_geom);

        KRATOS_ERROR_IF_NOT(p_brep)
            << "Geometry with id " << patch_id << " is not a BrepSurface\n";

        auto p_background =
            p_brep->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

        std::vector<double> knot_u, knot_v;
        p_background->SpansLocalSpace(knot_u, 0);
        p_background->SpansLocalSpace(knot_v, 1);

        KRATOS_ERROR_IF(knot_u.empty() || knot_v.empty())
            << "Empty knot vectors for patch " << patch_id << "\n";

        const double u_0  = knot_u.front();
        const double u_1  = knot_u.back();
        const double v_0 = knot_v.front();
        const double v_1 = knot_v.back();

        const IndexType number_pts = n_div + 1;

        pc.number_pts = number_pts;
        pc.u_0 = u_0;
        pc.v_0 = v_0;
        pc.delta_u  = (u_1  - u_0)  / static_cast<double>(n_div);
        pc.delta_v = (v_1 - v_0) / static_cast<double>(n_div);

        pc.points.reserve(number_pts * number_pts);

        CoordinatesArrayType local = ZeroVector(3);
        CoordinatesArrayType phys  = ZeroVector(3);

        IndexType id = 0;
        for (IndexType i = 0; i < number_pts; ++i) {
            local[0] = pc.u_0 + pc.delta_u * static_cast<double>(i);
            for (IndexType j = 0; j < number_pts; ++j) {
                local[1] = pc.v_0 + pc.delta_v * static_cast<double>(j);

                p_brep->GlobalCoordinates(phys, local);

                pc.points.push_back(PointTypePointer(new PointType(
                    id++, phys[0], phys[1], phys[2]
                )));
            }
        }

        pc.p_bins = std::make_unique<DynamicBins>(pc.points.begin(), pc.points.end());

        cache.emplace(patch_id, std::move(pc));
    }

    return cache;
}

// Finds IGA patches with sampled points near a specified physical-space location.
// Returns all input patches as a fallback when the radius search finds no candidate.
std::vector<IndexType> IgaMappingIntersectionUtilities::GetPatchesWithProbableProjection(
    const std::vector<IndexType>& patches_id,
    const PatchCacheMap& rCache,
    const CoordinatesArrayType& element_center,
    const double search_radius)
{
    std::vector<IndexType> out;
    out.reserve(patches_id.size());

    PointType search_point(0, element_center[0], element_center[1], element_center[2]);

    constexpr int number_of_results = 1;
    std::vector<PointTypePointer> results(number_of_results);
    std::vector<double> distances(number_of_results);

    for (const IndexType patch_id : patches_id) {
        auto it = rCache.find(patch_id);
        KRATOS_ERROR_IF(it == rCache.end())
            << "Patch id " << patch_id << " not found in patch cache.\n";

        const auto& pc = it->second;
        KRATOS_ERROR_IF_NOT(pc.p_bins)
            << "Bins not built for patch id " << patch_id << "\n";

        const int obtained = pc.p_bins->SearchInRadius(
            search_point,
            search_radius,
            results.begin(),
            distances.begin(),
            number_of_results);

        if (obtained > 0) out.push_back(patch_id);
    }

    return out.empty() ? patches_id : out;
}

// Creates surface quadrature-point coupling conditions from IGA-FEM coupling
// geometries. FEM triangles are projected, clipped, subdivided at knot spans,
// and integrated in parameter space.
void IgaMappingIntersectionUtilities::CreateIgaFEMQuadraturePointsOnSurface(
    ModelPart& rModelPartCoupling,
    bool origin_is_iga,
    const PatchCacheMap& rPatchCache,
    const double search_radius,
    const bool WriteTrianglesToFile)
{
    const ModelPart& rParentModelPart = rModelPartCoupling.GetParentModelPart();

    // Open the diagnostic output only when explicitly requested.
    std::ofstream parametric_triangles_output_file;
    std::ofstream physical_triangles_output_file;
    if (WriteTrianglesToFile) {
        const std::string parametric_output_file_name =
            "obtained_triangles_after_intersection_with_knot_lines_parameter_space.txt";
        const std::string physical_output_file_name =
            "obtained_triangles_after_intersection_with_knot_lines_physical_space.txt";

        parametric_triangles_output_file.open(parametric_output_file_name);
        physical_triangles_output_file.open(physical_output_file_name);
        KRATOS_ERROR_IF_NOT(parametric_triangles_output_file.is_open())
            << "Could not open triangle output file: "
            << parametric_output_file_name << "\n";
        KRATOS_ERROR_IF_NOT(physical_triangles_output_file.is_open())
            << "Could not open triangle output file: "
            << physical_output_file_name << "\n";

        parametric_triangles_output_file << std::setprecision(17)
            << "# Parametric triangle vertices: u v w\n"
            << "# Each triangle starts with: # brep_id <id>\n"
            << "# Every three non-empty lines define one triangle.\n";
        physical_triangles_output_file << std::setprecision(17)
            << "# Physical triangle vertices: x y z\n"
            << "# Every three non-empty lines define one triangle.\n";
    }

    // Iterate over the coupling geometries and create the quadrature points
    for (auto geometry_itr = rModelPartCoupling.GeometriesBegin();
        geometry_itr != rModelPartCoupling.GeometriesEnd();
        ++geometry_itr)
    {
        if (geometry_itr->NumberOfGeometryParts() < 2){ continue;}

        auto geom_master = geometry_itr->pGetGeometryPart(0); // IGA Surface
        auto geom_slave = geometry_itr->pGetGeometryPart(1);  // Finite element
        auto& r_geom_slave = *geom_slave;

        // Get the Brep surface representing the IGA patch (NURBS surface + boundary loops)
        auto geom_master_cast = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(geom_master);
        const auto p_master_nurbs_surface = dynamic_pointer_cast<
            NurbsSurfaceGeometry<3, PointerVector<Node>>>(
                geom_master_cast->pGetGeometryPart(
                    GeometryType::BACKGROUND_GEOMETRY_INDEX));
        KRATOS_ERROR_IF_NOT(p_master_nurbs_surface)
            << "The BRep background geometry is not a NURBS surface.\n";

        // This vector stores all the triangles resulting from the projection of a triangle to the surface patch (projection + subdivision)
        std::vector<std::vector<CoordinatesArrayType>> triangles_param_space;

        CoordinatesArrayType local_parameter = ZeroVector(3);
        CoordinatesArrayType node_coordinate_xyz;

        const IndexType n_nodes = r_geom_slave.size();

        std::vector<CoordinatesArrayType> successful_nodes_xyz;
        std::vector<CoordinatesArrayType> failed_nodes_xyz;
        std::vector<CoordinatesArrayType> points_to_triangulate;

        successful_nodes_xyz.reserve(n_nodes);
        failed_nodes_xyz.reserve(n_nodes);
        points_to_triangulate.reserve(n_nodes);

        for (IndexType i = 0; i < n_nodes; ++i){
            const auto p_point = r_geom_slave.pGetPoint(i);
            node_coordinate_xyz = p_point->GetInitialPosition();

            // Provide a good initial guess (reuses local_parameter)
            const bool have_initial_guess =
            IgaMappingIntersectionUtilities::FindInitialGuessNewtonRaphsonProjection(
                node_coordinate_xyz,        // slave XYZ
                geometry_itr->GetGeometryPart(0),                // Brep surface geometry
                rPatchCache,                // PatchCacheMap (built once)
                local_parameter,            // output (u,v)
                search_radius               // bins search radius
            );

            if (!have_initial_guess) {
            // Fallback if bins find nothing
                local_parameter = ZeroVector(3);
            }

            const bool projection_ok = geom_master->ProjectionPointGlobalToLocalSpace(
                node_coordinate_xyz, local_parameter, 1e-5) == 1;
            
            if (projection_ok)
            {
                successful_nodes_xyz.push_back(node_coordinate_xyz);
                points_to_triangulate.push_back(local_parameter);
            }
            else
            {
                failed_nodes_xyz.push_back(node_coordinate_xyz);
            }
        }

        const IndexType projection_is_successful_count =
            static_cast<IndexType>(successful_nodes_xyz.size());
        
        if (projection_is_successful_count == 0){
            if (!RecoverUnprojectedTriangleIntersection(
                    *geom_master, r_geom_slave, *p_master_nurbs_surface,
                    rPatchCache, search_radius, points_to_triangulate)) {
                continue;
            }
            MappingTriangulationUtilities::TriangulatePolygonFan(
                points_to_triangulate, triangles_param_space);
        }
        else if (projection_is_successful_count == 1){
            KRATOS_ERROR_IF(points_to_triangulate.size() != 1)
                << "Expected points_to_triangulate.size()==1, got "
                << points_to_triangulate.size() << "\n";

            KRATOS_ERROR_IF(failed_nodes_xyz.size() < 2)
                << "Expected at least 2 failed nodes, got "
                << failed_nodes_xyz.size() << "\n";

            // Besides the projected vertex and the two segment intersections,
            // the intersection polygon can contain corners of the rectangular
            // NURBS parameter domain.
            points_to_triangulate.reserve(7);

            bool all_found = true;
            for (IndexType i = 0; i < 2; ++i)
            {
                CoordinatesArrayType intersection_point_local_space = ZeroVector(3);

                const bool found =
                    FindTriangleSegmentSurfaceIntersectionWithBisection(
                        geometry_itr->GetGeometryPart(0),
                        successful_nodes_xyz[0],
                        failed_nodes_xyz[i],
                        points_to_triangulate[0],
                        intersection_point_local_space);

                if (!found) {
                    all_found = false;
                    break;
                }

                points_to_triangulate.push_back(intersection_point_local_space);
            }

            if (!all_found) {continue;}

            AddEnclosedDomainCorners(
                *p_master_nurbs_surface,
                r_geom_slave,
                points_to_triangulate);
            MappingTriangulationUtilities::TriangulatePolygonFan(
                points_to_triangulate, triangles_param_space);
        } else if (projection_is_successful_count == 2){
            KRATOS_ERROR_IF(points_to_triangulate.size() != 2)
                << "Expected 2 param points in points_to_triangulate, got "
                << points_to_triangulate.size() << "\n";

            KRATOS_ERROR_IF(failed_nodes_xyz.size() < 1)
                << "Expected at least 1 non-successful node.\n";

            // Two projected vertices, two segment intersections, and possibly
            // corners of the rectangular NURBS parameter domain.
            points_to_triangulate.reserve(8);

            bool all_found = true;
            for (IndexType i = 0; i < 2; ++i)
            {
                CoordinatesArrayType intersection_point_local_space = ZeroVector(3);

                const bool found =
                    IgaMappingIntersectionUtilities::FindTriangleSegmentSurfaceIntersectionWithBisection(
                    geometry_itr->GetGeometryPart(0),
                    successful_nodes_xyz[i],
                    failed_nodes_xyz[0],
                    points_to_triangulate[i],
                    intersection_point_local_space);

                if (!found) {
                    all_found = false;
                    break;
                }

                AddUniquePoint(
                    points_to_triangulate,
                    intersection_point_local_space);
            }

            if (!all_found) {continue;}

            AddEnclosedDomainCorners(
                *p_master_nurbs_surface,
                r_geom_slave,
                points_to_triangulate);
            MappingTriangulationUtilities::TriangulatePolygonFan(
                points_to_triangulate, triangles_param_space);
        } else if (projection_is_successful_count == 3)
        {
            SortVerticesCounterClockwise(points_to_triangulate);
            triangles_param_space.push_back(points_to_triangulate);
        } else {
            KRATOS_WARNING("IgaMappingIntersectionUtilities")
                << "The FEM-IGA projection resulted in a failure";
            continue;
        }

        // Clip every candidate triangle, irrespective of how many FEM nodes
        // were initially projected inside the trimmed surface.
        constexpr double factor = 1e-10;
        triangles_param_space =
            MappingTriangulationUtilities::ClipTrianglesWithTrimmingLoops(
            triangles_param_space,
            *geom_master_cast,
            factor);

        std::vector<MappingTriangulationUtilities::TriangleType> obtained_triangles_after_intersection_with_knot_lines;

        // Reuse quadrature containers outside loops
        IntegrationPointsArrayType integration_points_master(3);
        GeometriesArrayType quadrature_point_geometries_master(3);

        IntegrationPointsArrayType integration_points_slave(3);
        GeometriesArrayType quadrature_point_geometries_slave(3);

        CoordinatesArrayType master_quadrature_point_xyz = ZeroVector(3);
        CoordinatesArrayType slave_quadrature_point_local_space = ZeroVector(3);

        IntegrationInfo master_integration_info = geom_master->GetDefaultIntegrationInfo();

        for (IndexType tri_id = 0; tri_id < triangles_param_space.size(); ++tri_id)
        {
            obtained_triangles_after_intersection_with_knot_lines.clear();

            KRATOS_ERROR_IF(triangles_param_space[tri_id].size() != 3)
                << "Expected a triangle with 3 vertices, got "
                << triangles_param_space[tri_id].size() << ".\n";

            const MappingTriangulationUtilities::TriangleType triangle{{
                triangles_param_space[tri_id][0],
                triangles_param_space[tri_id][1],
                triangles_param_space[tri_id][2]}};

            MappingTriangulationUtilities::Triangulation(
                triangle, geom_master,
                obtained_triangles_after_intersection_with_knot_lines);

            for (IndexType j = 0; j < obtained_triangles_after_intersection_with_knot_lines.size(); ++j)
            {
                const auto& r_triangle = obtained_triangles_after_intersection_with_knot_lines[j];

                // Writing is optional and does not affect quadrature-point creation.
                if (WriteTrianglesToFile) {
                    WriteParametricTriangleToTextFile(
                        parametric_triangles_output_file,
                        r_triangle,
                        geom_master_cast->Id());
                    WritePhysicalTriangleToTextFile(
                        physical_triangles_output_file,
                        r_triangle,
                        *geom_master_cast);
                }

                const double xi_0  = r_triangle[0][0];
                const double xi_1  = r_triangle[1][0];
                const double xi_2  = r_triangle[2][0];
                const double eta_0 = r_triangle[0][1];
                const double eta_1 = r_triangle[1][1];
                const double eta_2 = r_triangle[2][1];

                auto master_it = integration_points_master.begin();
                IntegrationPointUtilities::IntegrationPointsTriangle2D(
                    master_it, 1, xi_0, xi_1, xi_2, eta_0, eta_1, eta_2);

                geom_master->CreateQuadraturePointGeometries(
                    quadrature_point_geometries_master, 3, integration_points_master, master_integration_info);

                integration_points_slave = integration_points_master;

                for (IndexType qp = 0; qp < quadrature_point_geometries_master.size(); ++qp)
                {
                    master_quadrature_point_xyz = quadrature_point_geometries_master[qp].Center();

                    geom_slave->PointLocalCoordinates(
                        slave_quadrature_point_local_space, master_quadrature_point_xyz);

                    integration_points_slave[qp].X() = slave_quadrature_point_local_space[0];
                    integration_points_slave[qp].Y() = slave_quadrature_point_local_space[1];
                }

                CreateQuadraturePointsUtility<NodeType>::Create(
                    r_geom_slave, quadrature_point_geometries_slave, integration_points_slave, 1);

                IndexType base_id =
                    (rParentModelPart.NumberOfConditions() == 0)
                        ? 1
                        : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;

                if (origin_is_iga)
                {
                    for (IndexType qp = 0; qp < quadrature_point_geometries_master.size(); ++qp)
                    {
                        rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                            base_id + qp,
                            Kratos::make_shared<CouplingGeometry<Node>>(
                                quadrature_point_geometries_master(qp),
                                quadrature_point_geometries_slave(qp))));
                    }
                }
                else
                {
                    for (IndexType qp = 0; qp < quadrature_point_geometries_master.size(); ++qp)
                    {
                        rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                            base_id + qp,
                            Kratos::make_shared<CouplingGeometry<Node>>(
                                quadrature_point_geometries_slave(qp),
                                quadrature_point_geometries_master(qp))));
                    }
                }
            }
        }
    }
}

// Finds a nearby cached patch sample and converts it into an initial parametric
// projection guess. Returns true when a cached point is found within the search radius.
bool IgaMappingIntersectionUtilities::FindInitialGuessNewtonRaphsonProjection(
    const CoordinatesArrayType& slave_xyz,
    const GeometryType& r_master_geometry,
    const PatchCacheMap& rPatchCache,
    CoordinatesArrayType& initial_guess,
    const double search_radius)
{
    const IndexType patch_id = r_master_geometry.Id();

    auto it = rPatchCache.find(patch_id);
    KRATOS_ERROR_IF(it == rPatchCache.end())
        << "Patch id " << patch_id << " not found in patch cache.\n";

    const PatchSearchCache& pc = it->second;
    KRATOS_ERROR_IF_NOT(pc.p_bins)
        << "Bins not built for patch id " << patch_id << "\n";

    PointType search_point(0, slave_xyz[0], slave_xyz[1], slave_xyz[2]);

    constexpr int number_of_results = 1;
    std::vector<PointTypePointer> results(number_of_results);
    std::vector<double> distances(number_of_results);

    const int obtained = pc.p_bins->SearchInRadius(
        search_point,
        search_radius,
        results.begin(),
        distances.begin(),
        number_of_results);

    if (obtained <= 0 || !results[0]) return false;

    const IndexType id = results[0]->Id();

    // Map id -> (i,j)
    const IndexType i = id / pc.number_pts;
    const IndexType j = id % pc.number_pts;

    initial_guess = ZeroVector(3);
    initial_guess[0] = pc.u_0  + pc.delta_u  * static_cast<double>(i);
    initial_guess[1] = pc.v_0 + pc.delta_v * static_cast<double>(j);

    return true;
}

// Checks whether projected points lie on a boundary of the NURBS parameter domain.
// Returns true for a boundary point or for points sharing the same boundary.
bool IgaMappingIntersectionUtilities::AreProjectionsOnParameterSpaceBoundary(
    const std::vector<CoordinatesArrayType>& r_points_to_triangulate,
    const NurbsSurfaceGeometry<3, PointerVector<Node>>& r_nurbs_surface)
{
    if (r_points_to_triangulate.empty()) {
        return false;
    }

    const NurbsInterval interval_u = r_nurbs_surface.DomainIntervalU();
    const NurbsInterval interval_v = r_nurbs_surface.DomainIntervalV(); 

    const double u_min = interval_u.MinParameter();
    const double u_max = interval_u.MaxParameter();
    const double v_min = interval_v.MinParameter();
    const double v_max = interval_v.MaxParameter();

    constexpr double eps = 1.0e-10;

    auto near = [](double a, double b) noexcept {
        return std::abs(a - b) < eps;
    };

    const auto& p0 = r_points_to_triangulate[0];

    // One point: on any boundary
    if (r_points_to_triangulate.size() == 1) {
        return near(p0[0], u_min) || near(p0[0], u_max) ||
               near(p0[1], v_min) || near(p0[1], v_max);
    }

    // Two or more points: check first two points on same boundary edge
    const auto& p1 = r_points_to_triangulate[1];

    const bool on_u_min = near(p0[0], u_min) && near(p1[0], u_min);
    const bool on_u_max = near(p0[0], u_max) && near(p1[0], u_max);
    const bool on_v_min = near(p0[1], v_min) && near(p1[1], v_min);
    const bool on_v_max = near(p0[1], v_max) && near(p1[1], v_max);

    return on_u_min || on_u_max || on_v_min || on_v_max;
}

// Locates the patch-boundary intersection of an inside-outside segment by
// bisection. Returns true after computing the best available local coordinates.
bool IgaMappingIntersectionUtilities::FindTriangleSegmentSurfaceIntersectionWithBisection(
    const GeometryType& r_geom_master,
    const CoordinatesArrayType& r_point_inside,
    const CoordinatesArrayType& r_point_outside,
    const CoordinatesArrayType& r_initial_guess,
    CoordinatesArrayType& r_intersection_point)
{
    CoordinatesArrayType a = r_point_inside;
    CoordinatesArrayType b = r_point_outside;
    CoordinatesArrayType inside_local_parameter = r_initial_guess;

    constexpr IndexType max_iters = 1000;
    constexpr double    x_tol     = 1e-8;  // tolerance in physical space (segment length)
    constexpr double    proj_tol  = 1e-6;  // projection tolerance

    for (IndexType it = 0; it < max_iters; ++it)
    {
        const CoordinatesArrayType mid = 0.5 * (a + b);
        CoordinatesArrayType trial_local_parameter = inside_local_parameter;

        // "Inside test" via projection success
        const bool projection_ok = r_geom_master.ProjectionPointGlobalToLocalSpace(
            mid, trial_local_parameter, proj_tol) == 1;

        if (!projection_ok) {
            // midpoint behaves as outside -> move outside endpoint inward
            b = mid;
        } else {
            // midpoint behaves as inside -> move inside endpoint outward
            a = mid;
            inside_local_parameter = trial_local_parameter;
        }

        if (norm_2(b - a) < x_tol) {
            r_intersection_point = inside_local_parameter;
            return true;
        }
    }

    r_intersection_point = inside_local_parameter;
    return false;
}

// Sorts parametric polygon vertices counter-clockwise around their centroid.
void IgaMappingIntersectionUtilities::SortVerticesCounterClockwise(
    std::vector<CoordinatesArrayType>& r_vertices)
{
    KRATOS_ERROR_IF(r_vertices.size() < 3)
        << "SortVerticesCounterClockwise expects at least 3 vertices, got "
        << r_vertices.size() << "\n";

    // --- centroid in parametric space ---
    double cx = 0.0;
    double cy = 0.0;

    for (const auto& v : r_vertices) {
        cx += v[0];
        cy += v[1];
    }
    cx /= static_cast<double>(r_vertices.size());
    cy /= static_cast<double>(r_vertices.size());

    // --- sort by angle around centroid ---
    std::sort(
        r_vertices.begin(),
        r_vertices.end(),
        [cx, cy](const CoordinatesArrayType& a, const CoordinatesArrayType& b)
        {
            const double ang_a = std::atan2(a[1] - cy, a[0] - cx);
            const double ang_b = std::atan2(b[1] - cy, b[0] - cx);
            return ang_a < ang_b;
        }
    );
}

namespace
{

// Writes one parametric triangle and its BRep patch ID to a text stream. Each
// vertex is written on a separate line as u, v, and w, followed by an empty line.
// rOutput: Stream receiving the triangle coordinates.
// rTriangle: Container holding exactly three parametric vertices.
// BrepId: ID of the BRep patch containing the triangle.
template<class TTriangleType>
void WriteParametricTriangleToTextFile(
    std::ostream& rOutput,
    const TTriangleType& rTriangle,
    const std::size_t BrepId)
{
    KRATOS_ERROR_IF(rTriangle.size() != 3)
        << "Expected a triangle with 3 vertices, got " << rTriangle.size() << ".\n";

    rOutput << "# brep_id " << BrepId << '\n';
    for (const auto& r_vertex : rTriangle) {
        rOutput << r_vertex[0] << ' ' << r_vertex[1] << ' ' << r_vertex[2] << '\n';
    }
    rOutput << '\n';
}

// Maps one parametric triangle to physical space and writes its vertices.
// rOutput: Stream receiving the physical triangle coordinates.
// rTriangle: Container holding exactly three parametric vertices.
// rBrepGeometry: BRep geometry used for the parameter-to-physical mapping.
template<class TTriangleType>
void WritePhysicalTriangleToTextFile(
    std::ostream& rOutput,
    const TTriangleType& rTriangle,
    const GeometryType& rBrepGeometry)
{
    KRATOS_ERROR_IF(rTriangle.size() != 3)
        << "Expected a triangle with 3 vertices, got "
        << rTriangle.size() << ".\n";

    CoordinatesArrayType physical_vertex = ZeroVector(3);
    for (const auto& r_parametric_vertex : rTriangle) {
        rBrepGeometry.GlobalCoordinates(
            physical_vertex, r_parametric_vertex);
        rOutput << physical_vertex[0] << ' '
                << physical_vertex[1] << ' '
                << physical_vertex[2] << '\n';
    }
    rOutput << '\n';
}

// Adds a parametric point unless an equivalent point is already present.
// rPoints: Point collection that can be extended by the function.
// rPoint: Candidate point to insert.
// Returns true when the point is inserted and false when it is a duplicate.
bool AddUniquePoint(
    std::vector<CoordinatesArrayType>& rPoints,
    const CoordinatesArrayType& rPoint)
{
    const bool already_present = std::any_of(
        rPoints.begin(), rPoints.end(),
        [&rPoint](const CoordinatesArrayType& rExistingPoint) {
            return norm_2(rExistingPoint - rPoint) <
                parametric_point_tolerance;
        });

    if (!already_present) {
        rPoints.push_back(rPoint);
    }
    return !already_present;
}

// Checks whether a physical point lies inside the initial FEM triangle. The test
// evaluates the barycentric coordinates using the initial positions of the three
// FEM nodes and includes points on the boundary within a numerical tolerance.
// rPoint: Physical point to test.
// rFEMGeometry: Three-node FEM triangle used for the containment test.
// Returns true when the point lies inside or on the triangle boundary.
bool IsInsideInitialFEMTriangle(
    const CoordinatesArrayType& rPoint,
    const GeometryType& rFEMGeometry)
{
    KRATOS_ERROR_IF(rFEMGeometry.size() != 3)
        << "Expected a three-node FEM triangle, got "
        << rFEMGeometry.size() << " nodes.\n";

    const CoordinatesArrayType a =
        rFEMGeometry.pGetPoint(0)->GetInitialPosition();
    const CoordinatesArrayType ab =
        rFEMGeometry.pGetPoint(1)->GetInitialPosition() - a;
    const CoordinatesArrayType ac =
        rFEMGeometry.pGetPoint(2)->GetInitialPosition() - a;
    const CoordinatesArrayType ap = rPoint - a;

    const double ab_ab = inner_prod(ab, ab);
    const double ab_ac = inner_prod(ab, ac);
    const double ac_ac = inner_prod(ac, ac);
    const double denominator = ab_ab * ac_ac - ab_ac * ab_ac;
    KRATOS_ERROR_IF(
        std::abs(denominator) <= std::numeric_limits<double>::epsilon())
        << "Cannot test a point against a degenerate FEM triangle.\n";

    const double coordinate_1 =
        (ac_ac * inner_prod(ap, ab) - ab_ac * inner_prod(ap, ac)) /
        denominator;
    const double coordinate_2 =
        (ab_ab * inner_prod(ap, ac) - ab_ac * inner_prod(ap, ab)) /
        denominator;
    constexpr double inside_tolerance = 1e-6;

    return coordinate_1 >= -inside_tolerance &&
           coordinate_2 >= -inside_tolerance &&
           coordinate_1 + coordinate_2 <= 1.0 + inside_tolerance;
}

// Adds NURBS parameter-domain corners enclosed by a FEM triangle. The four
// rectangular parameter-domain corners are mapped to physical space, tested
// against the initial FEM triangle, and added without duplicates.
// rNurbsSurface: NURBS surface providing the parameter domain and mapping.
// rFEMGeometry: FEM triangle tested for corner containment.
// rPoints: Parametric intersection polygon extended with enclosed corners.
void AddEnclosedDomainCorners(
    const NurbsSurfaceType& rNurbsSurface,
    const GeometryType& rFEMGeometry,
    std::vector<CoordinatesArrayType>& rPoints)
{
    const NurbsInterval interval_u = rNurbsSurface.DomainIntervalU();
    const NurbsInterval interval_v = rNurbsSurface.DomainIntervalV();
    const std::array<CoordinatesArrayType, 4> domain_corners{{
        CoordinatesArrayType{interval_u.MinParameter(), interval_v.MinParameter(), 0.0},
        CoordinatesArrayType{interval_u.MaxParameter(), interval_v.MinParameter(), 0.0},
        CoordinatesArrayType{interval_u.MaxParameter(), interval_v.MaxParameter(), 0.0},
        CoordinatesArrayType{interval_u.MinParameter(), interval_v.MaxParameter(), 0.0}}};

    for (const auto& r_corner_local : domain_corners) {
        CoordinatesArrayType corner_global = ZeroVector(3);
        rNurbsSurface.GlobalCoordinates(corner_global, r_corner_local);
        if (IsInsideInitialFEMTriangle(corner_global, rFEMGeometry)) {
            AddUniquePoint(rPoints, r_corner_local);
        }
    }
}

// Recovers an intersection when no FEM vertex projects onto the surface. The
// function adds enclosed NURBS-domain corners, samples every FEM edge to find
// projectable interior points, and uses bisection toward both edge endpoints to
// recover the parameter-space boundary intersections.
// rMasterGeometry: IGA geometry used for projection and bisection.
// rFEMGeometry: Three-node FEM triangle whose intersection is recovered.
// rNurbsSurface: NURBS surface used to add enclosed domain corners.
// rPatchCache: Spatial cache used to initialize surface projections.
// SearchRadius: Radius used when searching the projection cache.
// rPolygon: Recovered parametric polygon vertices.
// Returns true when at least three polygon vertices are recovered.
bool RecoverUnprojectedTriangleIntersection(
    const GeometryType& rMasterGeometry,
    const GeometryType& rFEMGeometry,
    const NurbsSurfaceType& rNurbsSurface,
    const PatchCacheMap& rPatchCache,
    const double SearchRadius,
    std::vector<CoordinatesArrayType>& rPolygon)
{
    AddEnclosedDomainCorners(rNurbsSurface, rFEMGeometry, rPolygon);

    constexpr std::size_t number_of_edge_samples = 16;
    for (std::size_t i = 0; i < 3; ++i) {
        const CoordinatesArrayType edge_point_1 =
            rFEMGeometry.pGetPoint(i)->GetInitialPosition();
        const CoordinatesArrayType edge_point_2 =
            rFEMGeometry.pGetPoint((i + 1) % 3)->GetInitialPosition();

        CoordinatesArrayType projected_edge_point = ZeroVector(3);
        CoordinatesArrayType projected_edge_local = ZeroVector(3);
        bool edge_crosses_domain = false;

        for (std::size_t sample = 1; sample < number_of_edge_samples; ++sample) {
            const double position = static_cast<double>(sample) /
                static_cast<double>(number_of_edge_samples);
            const CoordinatesArrayType edge_point =
                (1.0 - position) * edge_point_1 + position * edge_point_2;

            CoordinatesArrayType edge_local = ZeroVector(3);
            IgaMappingIntersectionUtilities::FindInitialGuessNewtonRaphsonProjection(
                edge_point, rMasterGeometry, rPatchCache,
                edge_local, SearchRadius);

            if (rMasterGeometry.ProjectionPointGlobalToLocalSpace(
                    edge_point, edge_local, 1e-5) == 1) {
                projected_edge_point = edge_point;
                projected_edge_local = edge_local;
                edge_crosses_domain = true;
                break;
            }
        }

        if (!edge_crosses_domain) {
            continue;
        }

        for (const auto& r_endpoint : {edge_point_1, edge_point_2}) {
            CoordinatesArrayType intersection_local = ZeroVector(3);
            if (IgaMappingIntersectionUtilities::
                    FindTriangleSegmentSurfaceIntersectionWithBisection(
                        rMasterGeometry, projected_edge_point, r_endpoint,
                        projected_edge_local, intersection_local)) {
                AddUniquePoint(rPolygon, intersection_local);
            }
        }
    }

    return rPolygon.size() >= 3;
}

} // unnamed namespace

} // namespace Kratos.
