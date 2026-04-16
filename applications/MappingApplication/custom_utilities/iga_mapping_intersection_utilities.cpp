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

// External includes

// Project includes
#include "iga_mapping_intersection_utilities.h"

// Data structures for doing spatial search 
#include "utilities/function_parser_utility.h"
#include "spatial_containers/bins_dynamic.h"

namespace Kratos
{

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


void IgaMappingIntersectionUtilities::CreateIgaFEMQuadraturePointsOnSurface(
    ModelPart& rModelPartCoupling,
    bool origin_is_iga,
    const PatchCacheMap& rPatchCache,
    const double search_radius)
{
    const ModelPart& rParentModelPart = rModelPartCoupling.GetParentModelPart();

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

        // Get the NURBS Surface defining the Brep surface 
        GeometryPointerType master_nurbs_surface = geom_master_cast->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
        auto master_nurbs_surface_cast = dynamic_pointer_cast<NurbsSurfaceGeometry<3, PointerVector<Node>>>(master_nurbs_surface);

        // This vector stores all the triangles resulting from the projection of a triangle to the surface patch (projection + subdivision)
        std::vector<std::vector<CoordinatesArrayType>> triangles_param_space;

        CoordinatesArrayType local_parameter = ZeroVector(3);
        CoordinatesArrayType node_coordinate_xyz;

        const IndexType n_nodes = r_geom_slave.size();

        // These nodes refer to the FEM nodes which are being projected to the IGA patch
        std::vector<CoordinatesArrayType> successful_nodes_xyz;
        std::vector<CoordinatesArrayType> failed_nodes_xyz;
        std::vector<CoordinatesArrayType> points_to_triangulate;

        // Reserve (worst case: all nodes succeed or fail)
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
        
        // If no node is projected inside the patch, continue to the next iteration. Otherwise, find the intersection between the triangle sides and the surface
        if (projection_is_successful_count == 0){
            continue;
        }
        else if (projection_is_successful_count == 1){  // If just one node is projected inside the surface, create one triangle  
            
            KRATOS_ERROR_IF(points_to_triangulate.size() != 1)
                << "Expected points_to_triangulate.size()==1, got "
                << points_to_triangulate.size() << "\n";

            KRATOS_ERROR_IF(failed_nodes_xyz.size() < 2)
                << "Expected at least 2 failed nodes, got "
                << failed_nodes_xyz.size() << "\n";

            // Check if the projection is on the boundary of the parameter space and, if so, continue to the next element
            if (AreProjectionsOnParameterSpaceBoundary(points_to_triangulate, *master_nurbs_surface_cast) == true) continue;

            // Vector for storing the intersection of the triangle side with the boundaries of the patch
            std::vector<CoordinatesArrayType> triangle_segment_intersection_with_surface_patch_local_parameter;
            triangle_segment_intersection_with_surface_patch_local_parameter.reserve(2);

            points_to_triangulate.reserve(3); // 1 existing + 2 intersections

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

                triangle_segment_intersection_with_surface_patch_local_parameter.push_back(intersection_point_local_space);
                points_to_triangulate.push_back(intersection_point_local_space);
            }

            if (!all_found) {continue;}

            // If the intersection between the triangle and the patch surface is at just one point, ignore the element and go to the next triangle
            if (norm_2(triangle_segment_intersection_with_surface_patch_local_parameter[0] - triangle_segment_intersection_with_surface_patch_local_parameter[1]) < 1e-4)
            {
                continue;
            }

            KRATOS_ERROR_IF(points_to_triangulate.size() != 3)
                << "Expected 3 points to form a triangle, got "
                << points_to_triangulate.size() << "\n";

            // Sort the vertices counter-clockwise before doing triangulation
            SortVerticesCounterClockwise(points_to_triangulate);

            // Form one triangle with the one point projected inside the patch and the two intersections
            triangles_param_space.push_back(points_to_triangulate);
        } else if (projection_is_successful_count == 2){  // If just one node is projected inside the surface, create 2 triangles

            // Skip if projections lie on parameter-space boundary
            if (AreProjectionsOnParameterSpaceBoundary(points_to_triangulate, *master_nurbs_surface_cast) == true) continue;

            KRATOS_ERROR_IF(points_to_triangulate.size() != 2)
                << "Expected 2 param points in points_to_triangulate, got "
                << points_to_triangulate.size() << "\n";

            KRATOS_ERROR_IF(failed_nodes_xyz.size() < 1)
                << "Expected at least 1 non-successful node.\n";

            // We will add up to 2 intersection points
            points_to_triangulate.reserve(4);

            constexpr double tol_dup = 1e-4;

            auto IsDuplicateParamPoint = [&](const CoordinatesArrayType& p) {
                for (const auto& q : points_to_triangulate) {
                    if (norm_2(p - q) < tol_dup) return true;
                }
                return false;
            };

            for (IndexType i = 0; i < 2; ++i)
            {
                CoordinatesArrayType intersection_point_local_space = ZeroVector(3);

                // Use your current intersection routine (void). If you have the new bool routine, use it instead.
                IgaMappingIntersectionUtilities::FindTriangleSegmentSurfaceIntersectionWithBisection(
                    geometry_itr->GetGeometryPart(0),
                    successful_nodes_xyz[i],              // inside the patch
                    failed_nodes_xyz[0],          // outside the patch
                    points_to_triangulate[i],                                       // better initial guess: the corresponding inside param
                    intersection_point_local_space);

                if (!IsDuplicateParamPoint(intersection_point_local_space)) {
                    points_to_triangulate.push_back(intersection_point_local_space);
                }
            }

            SortVerticesCounterClockwise(points_to_triangulate);

            if (points_to_triangulate.size() == 3)
            {
                // Create one triangle
                triangles_param_space.push_back(points_to_triangulate);
            }
            else if (points_to_triangulate.size() == 4)
            {
                // Create two triangles: (0,1,2) and (0,2,3)
                triangles_param_space.push_back(
                    {points_to_triangulate[0], points_to_triangulate[1], points_to_triangulate[2]}
                );
                triangles_param_space.push_back(
                    {points_to_triangulate[0], points_to_triangulate[2], points_to_triangulate[3]}
                );
            }
        } else if (projection_is_successful_count == 3)
        { // If the 3 nodes are projected inside the untrimmed surface, check if the triangles intersect the trimming curves and , if so, create new triangles

            // Get the trimming curves outer loop
            auto outer_loop_array = geom_master_cast->GetOuterLoops();

            KRATOS_ERROR_IF(geom_master_cast->GetInnerLoops().size() > 0) << "The current patch has inner loops and this is not yet supported for the IGA-FEM Surface Mortar Mapper" << std::endl;

            Clipper2Lib::Paths64 all_loops(1), triangle(1), solution;
            const double factor = 1e-10;

            Clipper2Lib::Point64 int_point;
            int_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
            int_point.y = static_cast<cInt>(std::numeric_limits<int>::min());

            // Reuse tessellation object 
            CurveTessellation<PointerVector<Node>> curve_tessellation;

            // Reserve some space for the curve tessellation
            all_loops[0].reserve(all_loops[0].size() + 1000);

            for (IndexType j = 0; j < outer_loop_array[0].size(); ++j)
            {
                const auto& geometry_outer = *(outer_loop_array[0][j]);

                curve_tessellation.Tessellate(
                    geometry_outer, 0.001, 1, true);

                // Avoid copy
                const auto& tessellation = curve_tessellation.GetTessellation();

                for (IndexType u = 0; u < tessellation.size(); ++u)
                {
                    // Bind once
                    const auto& uv = std::get<1>(tessellation[u]);

                    auto new_int_point = BrepTrimmingUtilities<false>::ToIntPoint(
                        uv[0], uv[1], factor);

                    // We need to define this tolerance for the clipping process
                    new_int_point.x += 100000.0;
                    new_int_point.y += 100000.0;

                    // Fast duplicate rejection
                    if (new_int_point.x != int_point.x ||
                        new_int_point.y != int_point.y)
                    {
                        all_loops[0].push_back(new_int_point);
                        int_point = new_int_point;
                    }
                }
            }

            triangle[0].reserve(triangle[0].size() + 3);

            for (IndexType u = 0; u < 3; ++u){
                const auto new_int_point =
                    BrepTrimmingUtilities<false>::ToIntPoint(
                            points_to_triangulate[u][0],
                            points_to_triangulate[u][1],
                            factor);

                if (new_int_point.x != int_point.x ||
                    new_int_point.y != int_point.y)
                {
                    triangle[0].push_back(new_int_point);
                    int_point = new_int_point;
                }
            }

            // Intersect the triangle with the polygon resulting from the tessellation of the outer loop
            solution = Clipper2Lib::Intersect(all_loops, triangle, Clipper2Lib::FillRule::NonZero);

            double clip_area = 0.0;
            if (solution.size() > 0){
                clip_area = std::abs(Clipper2Lib::Area(solution[0]));
            }

            std::vector<Matrix> triangles;
            std::vector<CoordinatesArrayType> sorted_points_to_triangulate;
            sorted_points_to_triangulate.resize(3);

            if (!solution.empty()){
                triangles.clear();

                BrepTrimmingUtilities<false>::Triangulate_OPT(solution[0], triangles, factor, clip_area);

                triangles_param_space.reserve(
                    triangles_param_space.size() + triangles.size());

                for (IndexType u = 0; u < triangles.size(); ++u){
                    // Fill local triangle points
                    for (IndexType v = 0; v < 3; ++v){
                        points_to_triangulate[v][0] = triangles[u](v, 0);
                        points_to_triangulate[v][1] = triangles[u](v, 1);
                        points_to_triangulate[v][2] = 0.0;
                    }

                    // Reuse preallocated output container
                    SortVerticesCounterClockwise(points_to_triangulate);

                    // Avoid copy if possible
                    triangles_param_space.push_back(points_to_triangulate);
                }
                } else {    
                    KRATOS_WARNING("The intersection resulted in an empty polygon, skipping the triangle.");
                }

        } else {
            KRATOS_WARNING("IgaMappingIntersectionUtilities")
                << "The FEM-IGA projection resulted in a failure";
            continue;
        }

        std::vector<std::vector<CoordinatesArrayType>> obtained_triangles_after_intersection_with_knot_lines;

        // Reuse quadrature containers outside loops
        IntegrationPointsArrayType integration_points_master(3);
        GeometriesArrayType quadrature_point_geometries_master(3);

        IntegrationPointsArrayType integration_points_slave(3);
        GeometriesArrayType quadrature_point_geometries_slave(3);

        CoordinatesArrayType master_quadrature_point_physical_space = ZeroVector(3);
        CoordinatesArrayType slave_quadrature_point_local_space = ZeroVector(3);

        IntegrationInfo master_integration_info = geom_master->GetDefaultIntegrationInfo();

        for (IndexType tri_id = 0; tri_id < triangles_param_space.size(); ++tri_id)
        {
            obtained_triangles_after_intersection_with_knot_lines.clear();

            // Fallback: no triangulation => process original triangle
            if (obtained_triangles_after_intersection_with_knot_lines.empty()) {
                obtained_triangles_after_intersection_with_knot_lines.push_back(triangles_param_space[tri_id]);
            }

            for (IndexType j = 0; j < obtained_triangles_after_intersection_with_knot_lines.size(); ++j)
            {
                const auto& r_triangle = obtained_triangles_after_intersection_with_knot_lines[j];

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
                    master_quadrature_point_physical_space = quadrature_point_geometries_master[qp].Center();

                    geom_slave->PointLocalCoordinates(
                        slave_quadrature_point_local_space, master_quadrature_point_physical_space);

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

void IgaMappingIntersectionUtilities::CreateIgaIgaQuadraturePointsCoupling2DGeometries3D(
        ModelPart &rModelPartCoupling)
{
    const ModelPart &rParentModelPart = rModelPartCoupling.GetParentModelPart();

    // Iterate over the coupling geometries and create the quadrature points
    for (auto geometry_itr = rModelPartCoupling.GeometriesBegin();
            geometry_itr != rModelPartCoupling.GeometriesEnd();
            ++geometry_itr)
    {
        // Don´t enter the for loop if the geometry is a brep curve on surface or a brep surface (remember that in order to get the brep surface we have to include in the coupling model part the geometries of origin and destination)
        if (geometry_itr->GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Coupling_Geometry)
        {
            continue;
        }

        // Extract the master and slave parts from the coupling geometry
        auto geom_master = geometry_itr->pGetGeometryPart(0); // IGA patch
        auto geom_slave = geometry_itr->pGetGeometryPart(1);  // IGA patch

        //     auto &r_geom_slave = geometry_itr->GetGeometryPart(1);

        auto geom_master_cast = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(geom_master);
        auto geom_slave_cast = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, false, PointerVector<Point>>>(geom_slave);

        // Get the knot vectors in xi and eta direction for the destination domain
        std::vector<double> knot_vector_xi_slave_domain, knot_vector_eta_slave_domain;
        geom_slave->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->SpansLocalSpace(knot_vector_xi_slave_domain, 0);
        geom_slave->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->SpansLocalSpace(knot_vector_eta_slave_domain, 1);

        // Get the trimming curves outer loop
        auto outer_loop_array = geom_slave_cast->GetOuterLoops();

        // Construct the clipper library objects
        Clipper2Lib::Paths64 all_loops(1), rectangle(1), solution;
        const double factor = 1e-10;

        Clipper2Lib::Point64 int_point;
        int_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
        int_point.y = static_cast<cInt>(std::numeric_limits<int>::min());

        // Loop over the trimming curves, tesselate them and store the discretized geometry in the all_loops object
        for (IndexType j = 0; j < outer_loop_array[0].size(); ++j)
        {
            CurveTessellation<PointerVector<Node>> curve_tesselation;
            auto geometry_outer = *(outer_loop_array[0][j].get());
            curve_tesselation.Tessellate(
                geometry_outer, 0.001, 1, true);
            auto tesselation = curve_tesselation.GetTessellation();
            for (IndexType u = 0; u < tesselation.size(); ++u)
            {
                auto new_int_point = BrepTrimmingUtilities<false>::ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);
                if (!(int_point.x == new_int_point.x && int_point.y == new_int_point.y))
                {
                    all_loops[0].push_back(new_int_point);
                    int_point.x = new_int_point.x;
                    int_point.y = new_int_point.y;
                }
            }
        }

        // Loop over the "elements" in the destination domain
        std::vector<std::vector<CoordinatesArrayType>> all_new_rectangles_master_parameter_space;
        for (IndexType i = 0; i < knot_vector_xi_slave_domain.size() - 1; ++i)
        {
            for (IndexType j = 0; j < knot_vector_eta_slave_domain.size() - 1; ++j)
            {
                Clipper2Lib::Rect64 rectangle = Clipper2Lib::Rect64(
                    static_cast<cInt>(knot_vector_xi_slave_domain[i] / factor), static_cast<cInt>(knot_vector_eta_slave_domain[j] / factor),
                    static_cast<cInt>(knot_vector_xi_slave_domain[i + 1] / factor), static_cast<cInt>(knot_vector_eta_slave_domain[j + 1] / factor));

                // Intersect the rectangle with the polygon resulting from the tessellation of the outer loop
                solution = Clipper2Lib::RectClip(rectangle, all_loops);

                    // Area of the original rectangle
                    const double span_area = std::abs(Clipper2Lib::Area(rectangle.AsPath()));
                    double clip_area = 0.0;

                if (solution.size() > 0){
                    clip_area = std::abs(Clipper2Lib::Area(solution[0]));
                    for (IndexType k = 1; k < solution.size(); ++k){
                        clip_area -= std::abs(Clipper2Lib::Area(solution[k]));
                    }
                } else {
                    continue; // If the solution size is equal to 0, continue to the next iteration
                }

                // Check if the element is untrimmed. If it is untrimmed, map the rectangle to the destination physical space and then project to the origin parameter space
                if (std::abs(clip_area - span_area) < 1000){
                    // Create a vector with the coordinates of the rectangle in the destination physical space
                    std::vector<CoordinatesArrayType> knot_span_coordinates_physical_space, knot_span_coordinates_master_parameter_space;

                    // Map the rectangle vertices from the parameter space to the physical space
                    for (IndexType k = i; k < i + 2; k++){
                        for (IndexType h = j; h < j + 2; h++){
                            CoordinatesArrayType vertex_coordinate_physical_space, vertex_coordinate_slave_parameter_space{knot_vector_xi_slave_domain[k], knot_vector_eta_slave_domain[h], 0.0};

                            geom_slave->GlobalCoordinates(vertex_coordinate_physical_space, vertex_coordinate_slave_parameter_space);

                            knot_span_coordinates_physical_space.push_back(vertex_coordinate_physical_space);
                        }
                    }

                    // Project the rectangle to the origin domain parameter space
                    for (IndexType k = 0; k < knot_span_coordinates_physical_space.size(); k++){
                        CoordinatesArrayType vertex_coordinate_master_parameter_space;

                        geom_master->ProjectionPointGlobalToLocalSpace(knot_span_coordinates_physical_space[k], vertex_coordinate_master_parameter_space, 1e-6);

                        knot_span_coordinates_master_parameter_space.push_back(vertex_coordinate_master_parameter_space);
                    }

                    // Check if the projected rectangles are cut by knot lines in the origin parameter space and subdivide them in new rectangles
                    std::vector<std::vector<CoordinatesArrayType>> new_rectangles_master_parameter_space;
                    auto reordered_knot_span_coordinates_master_parameter_space = IgaMappingIntersectionUtilities::OrderRectangleCoordinates(knot_span_coordinates_master_parameter_space); 
                    IgaMappingIntersectionUtilities::SubdivideRectangleWithMasterPatchKnotLines(reordered_knot_span_coordinates_master_parameter_space, geom_master, new_rectangles_master_parameter_space);
                    all_new_rectangles_master_parameter_space.insert(
                        all_new_rectangles_master_parameter_space.end(),
                        new_rectangles_master_parameter_space.begin(),
                        new_rectangles_master_parameter_space.end());

                    // Iterate over the newly created rectangles in the origin parameter space
                    for (IndexType k = 0; k < new_rectangles_master_parameter_space.size(); k++){
                        // Create integration points in the master parameter space
                        IntegrationInfo master_integration_info = geom_master->GetDefaultIntegrationInfo();
                        IntegrationInfo slave_integration_info = geom_slave->GetDefaultIntegrationInfo();
                        const IndexType number_of_integration_points_master = master_integration_info.GetNumberOfIntegrationPointsPerSpan(0) * master_integration_info.GetNumberOfIntegrationPointsPerSpan(1);
                        const IndexType number_of_integration_points_slave = master_integration_info.GetNumberOfIntegrationPointsPerSpan(0) * master_integration_info.GetNumberOfIntegrationPointsPerSpan(1);
                        IndexType number_of_integration_points;

                        // Define the number of integration points to place in the origin domain as the biggest number between the origin and destination domains
                        number_of_integration_points = (number_of_integration_points_master >= number_of_integration_points_slave)
                            ? number_of_integration_points_master
                            : number_of_integration_points_slave;

                        // Create integration points and quadrature points containers for the master
                        IntegrationPointsArrayType integration_points_master(number_of_integration_points);
                        typename IntegrationPointsArrayType::iterator integration_point_master_iterator = integration_points_master.begin();
                        GeometriesArrayType quadrature_point_geometries_master(number_of_integration_points);

                        // Create integration and quadrature points containers for the slave
                        IntegrationPointsArrayType integration_points_slave(number_of_integration_points);
                        GeometriesArrayType quadrature_point_geometries_slave(number_of_integration_points);

                        IntegrationPointUtilities::IntegrationPoints2D(
                            integration_point_master_iterator,
                            master_integration_info.GetNumberOfIntegrationPointsPerSpan(0), master_integration_info.GetNumberOfIntegrationPointsPerSpan(1),
                            new_rectangles_master_parameter_space[k][0][0], new_rectangles_master_parameter_space[k][1][0],
                            new_rectangles_master_parameter_space[k][0][1], new_rectangles_master_parameter_space[k][2][1]);

                        // Create master quadrature points
                        geom_master->CreateQuadraturePointGeometries(quadrature_point_geometries_master, 3, integration_points_master, master_integration_info);

                        // Get the position of the master quadrature points in its physical space
                        CoordinatesArrayType master_quadrature_point_physical_space = ZeroVector(3);
                        CoordinatesArrayType slave_quadrature_point_local_space = ZeroVector(3);
                        integration_points_slave = integration_points_master;

                        for (IndexType k = 0; k < quadrature_point_geometries_master.size(); k++){

                            // Get the position of the master i_th quadrature point in its physical space
                            master_quadrature_point_physical_space = quadrature_point_geometries_master[k].Center();

                            // Obtain the local coordinate of the i_th master quadrature point in the slave local space
                            geom_slave->ProjectionPointGlobalToLocalSpace(master_quadrature_point_physical_space, slave_quadrature_point_local_space, 1e-6);

                            // Modify the slave integration points container
                            integration_points_slave[k].X() = slave_quadrature_point_local_space[0];
                            integration_points_slave[k].Y() = slave_quadrature_point_local_space[1];
                        }

                        // Create slave quadrature points
                        geom_slave->CreateQuadraturePointGeometries(quadrature_point_geometries_slave, 3, integration_points_slave, slave_integration_info);

                        // Create new conditions in the coupling model part whose geometry is CouplingGeometry
                        const IndexType id = (rParentModelPart.NumberOfConditions() == 0) ? 1 : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;

                        for (IndexType i = 0; i < quadrature_point_geometries_master.size(); ++i){
                            rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                                id + i, Kratos::make_shared<CouplingGeometry<Node>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i))));
                        }
        //             }
        //             else
        //             {
        //                 std::vector<std::vector<CoordinatesArrayType>> triangles_destination_parameter_space; // Vector of triangles coming from the triangulation process (="triangles", but of different kind)
        //                 std::vector<Matrix> triangles;

        //                 // Do triangulation
        //                 BrepTrimmingUtilities::Triangulate_OPT(solution[0], triangles, factor);

        //                 for (IndexType u = 0; u < triangles.size(); ++u)
        //                 {
        //                     std::vector<CoordinatesArrayType> triangle; // One specific triangle coming from the triangulation process
        //                     // Extract the triangles from the Triangulation_OPT method
        //                     for (IndexType vertex_index = 0; vertex_index < 3; vertex_index++)
        //                     {
        //                         double coordinate_xi = triangles[u](vertex_index, 0);
        //                         double coordinate_eta = triangles[u](vertex_index, 1);
        //                         CoordinatesArrayType node_coordinate{coordinate_xi, coordinate_eta, 0.0};
        //                         triangle.push_back(node_coordinate);
        //                     }

        //                     // Order the triangles vertices counterclockwise
        //                     std::vector<CoordinatesArrayType> sorted_triangle;
        //                     sortVerticesCounterClockwise(triangle, sorted_triangle);

        //                     triangles_destination_parameter_space.push_back(sorted_triangle);
        //                 }

        //                 for (IndexType k = 0; k < triangles_destination_parameter_space.size(); k++)
        //                 {

        //                     // Create a vector with the coordinates of the triangle in the destination physical space
        //                     std::vector<CoordinatesArrayType> triangle_coordinates_xyz, triangle_coordinates_master_parameter_space;

        //                     // Map each triangle node from parameter to physical space
        //                     for (IndexType h = 0; h < 3; h++)
        //                     {
        //                         CoordinatesArrayType vertex_coordinate_physical_space;
        //                         geom_slave->GlobalCoordinates(vertex_coordinate_physical_space, triangles_destination_parameter_space[k][h]);
        //                         triangle_coordinates_xyz.push_back(vertex_coordinate_physical_space);
        //                     }

        //                     // Project the triangle to the origin domain parameter space
        //                     for (IndexType k = 0; k < triangle_coordinates_xyz.size(); k++)
        //                     {
        //                         CoordinatesArrayType vertex_coordinate_master_parameter_space;

        //                         geom_master->ProjectionPointGlobalToLocalSpace(triangle_coordinates_xyz[k], vertex_coordinate_master_parameter_space, 1e-6);

        //                         triangle_coordinates_master_parameter_space.push_back(vertex_coordinate_master_parameter_space);
        //                     }

        //                     // Triangulation step: Intersect the triangles in the origin parameter space with the knot lines and subdivide them
        //                     std::vector<std::vector<CoordinatesArrayType>> obtained_triangles; // The outer vector contains the newly found triangles

        //                     // The triangulation process will create new triangles which are not intersected by the knot lines
        //                     IgaMappingIntersectionUtilities::Triangulation(triangle_coordinates_master_parameter_space, geom_master, obtained_triangles);

        //                     for (IndexType j = 0; j < obtained_triangles.size(); j++)
        //                     {
        //                         // Create integration points container for the master
        //                         IntegrationPointsArrayType integration_points_master(3);
        //                         typename IntegrationPointsArrayType::iterator integration_point_master_iterator = integration_points_master.begin();
        //                         GeometriesArrayType quadrature_point_geometries_master(3);

        //                         // Create integration points container for the slave
        //                         IntegrationPointsArrayType integration_points_slave(3);
        //                         GeometriesArrayType quadrature_point_geometries_slave(3);

        //                         // Coordinates of the triangle vertex in the master parameter space
        //                         double xi_0 = obtained_triangles[j][0][0];
        //                         double xi_1 = obtained_triangles[j][1][0];
        //                         double xi_2 = obtained_triangles[j][2][0];
        //                         double eta_0 = obtained_triangles[j][0][1];
        //                         double eta_1 = obtained_triangles[j][1][1];
        //                         double eta_2 = obtained_triangles[j][2][1];

        //                         // Creating the integration points
        //                         IntegrationPointUtilities::IntegrationPointsTriangle2D(integration_point_master_iterator, 1, xi_0, xi_1, xi_2, eta_0, eta_1, eta_2);

        //                         // Create master quadrature points
        //                         IntegrationInfo master_integration_info = geom_master->GetDefaultIntegrationInfo();
        //                         geom_master->CreateQuadraturePointGeometries(quadrature_point_geometries_master, 3, integration_points_master, master_integration_info);

        //                         // Get the position of the master quadrature points in its physical space
        //                         CoordinatesArrayType master_quadrature_point_physical_space = ZeroVector(3);
        //                         CoordinatesArrayType slave_quadrature_point_local_space = ZeroVector(3);
        //                         integration_points_slave = integration_points_master;

        //                         for (IndexType k = 0; k < quadrature_point_geometries_master.size(); k++)
        //                         {
        //                             // Get the position of the master i_th quadrature point in its physical space

        //                             master_quadrature_point_physical_space = quadrature_point_geometries_master[k].Center();

        //                             // Obtain the local coordinate of the i_th master quadrature point in the slave local space
        //                             geom_slave->ProjectionPointGlobalToLocalSpace(master_quadrature_point_physical_space, slave_quadrature_point_local_space, 1e-6);

        //                             // Modify the slave integration points container
        //                             integration_points_slave[k].X() = slave_quadrature_point_local_space[0];
        //                             integration_points_slave[k].Y() = slave_quadrature_point_local_space[1];
        //                         }

        //                         // Create slave quadrature points
        //                         IntegrationInfo slave_integration_info = geom_slave->GetDefaultIntegrationInfo();
        //                         geom_slave->CreateQuadraturePointGeometries(quadrature_point_geometries_slave, 3, integration_points_slave, slave_integration_info);

        //                         // add the quadrature point geometry conditions to the result model part
        //                         const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
        //                                                  ? 1
        //                                                  : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;

        //                         for (IndexType i = 0; i < quadrature_point_geometries_master.size(); ++i)
        //                         {
        //                             rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
        //                                 id + i, Kratos::make_shared<CouplingGeometry<Node>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i))));
        //                         }
        //                     }
                        }
                    }
                }
            }
            IgaMappingIntersectionUtilities::WriteRectanglesToPythonPlot(
                        all_new_rectangles_master_parameter_space,
                        "rectangles_plot.py");
        }
}


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

    constexpr double eps = 1.0e-2;

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

bool IgaMappingIntersectionUtilities::FindTriangleSegmentSurfaceIntersectionWithBisection(
    const GeometryType& r_geom_master,
    const CoordinatesArrayType& r_point_inside,
    const CoordinatesArrayType& r_point_outside,
    const CoordinatesArrayType& r_initial_guess,
    CoordinatesArrayType& r_intersection_point)
{
    CoordinatesArrayType a = r_point_inside;
    CoordinatesArrayType b = r_point_outside;

    // This will be updated by ProjectionPointGlobalToLocalSpace
    CoordinatesArrayType local_param = r_initial_guess;

    constexpr IndexType max_iters = 1000;
    constexpr double    x_tol     = 1e-8;  // tolerance in physical space (segment length)
    constexpr double    proj_tol  = 1e-6;  // projection tolerance (your original)

    CoordinatesArrayType mid = ZeroVector(3);

    for (IndexType it = 0; it < max_iters; ++it)
    {
        mid = 0.5 * (a + b);

        // Stop when segment is small enough (more meaningful than consecutive-mid distance)
        const double seg_len = norm_2(b - a);
        if (seg_len < x_tol) {
            r_intersection_point = local_param;
            return true;
        }

        // "Inside test" via projection success
        const IndexType projection_ok = r_geom_master.ProjectionPointGlobalToLocalSpace(
            mid, local_param, proj_tol);

        if (projection_ok == 0) {
            // midpoint behaves as outside -> move outside endpoint inward
            b = mid;
        } else {
            // midpoint behaves as inside -> move inside endpoint outward
            a = mid;
        }
    }

    // Best effort output after max iters
    r_intersection_point = local_param;
    return true; // or false if you want to signal "hit max iters"
}

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

void IgaMappingIntersectionUtilities::SubdivideRectangleWithMasterPatchKnotLines(
        const std::vector<CoordinatesArrayType>& r_original_knot_span_coordinates,
        GeometryPointerType master_geometry,
        std::vector<std::vector<CoordinatesArrayType>> &r_new_rectangles_master_parameter_space)
{
    r_new_rectangles_master_parameter_space.push_back(r_original_knot_span_coordinates);

    // Get the knot vectors of the master patch in the xi and eta directions
    std::vector<double> knot_vector_xi;
    std::vector<double> knot_vector_eta;

    master_geometry->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->SpansLocalSpace(knot_vector_xi, 0);
    master_geometry->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->SpansLocalSpace(knot_vector_eta, 1);

    // Check if both knot vectors size is 2 and if true, return the original rectangle coordinates (no new rectangles will be found)
    if (knot_vector_xi.size() == 2 && knot_vector_eta.size() == 2){ return; }

    // Obtain the inner knots for the xi and eta direction
    std::vector<double> inner_knots_xi, inner_knots_eta;

    if (knot_vector_xi.size() > 2){ inner_knots_xi = std::vector<double>(knot_vector_xi.begin() + 1, knot_vector_xi.end() - 1); }

    if (knot_vector_eta.size() > 2){ inner_knots_eta = std::vector<double>(knot_vector_eta.begin() + 1, knot_vector_eta.end() - 1); }

    for (IndexType k = 0; k < inner_knots_xi.size(); k++){

        std::vector<std::vector<CoordinatesArrayType>> intermediate_rectangles_xi;

        for (IndexType j = 0; j < r_new_rectangles_master_parameter_space.size(); j++){
            std::vector<std::vector<CoordinatesArrayType>> new_temp_rectangles_xi;
            SplitRectangle(new_temp_rectangles_xi, r_new_rectangles_master_parameter_space[j], inner_knots_xi[k], true);
            for (IndexType i = 0; i < new_temp_rectangles_xi.size(); i++){
                intermediate_rectangles_xi.push_back(new_temp_rectangles_xi[i]);
            }
        }

        r_new_rectangles_master_parameter_space = intermediate_rectangles_xi;
    }

    for (IndexType k = 0; k < inner_knots_eta.size(); k++){

        std::vector<std::vector<CoordinatesArrayType>> intermediate_rectangles_eta;

        for (IndexType j = 0; j < r_new_rectangles_master_parameter_space.size(); j++){
            std::vector<std::vector<CoordinatesArrayType>> new_temp_rectangles_eta;
            SplitRectangle(new_temp_rectangles_eta, r_new_rectangles_master_parameter_space[j], inner_knots_eta[k], false);
            for (IndexType i = 0; i < new_temp_rectangles_eta.size(); i++){
                intermediate_rectangles_eta.push_back(new_temp_rectangles_eta[i]);
            }
        }

        r_new_rectangles_master_parameter_space = intermediate_rectangles_eta;
    }
}

void IgaMappingIntersectionUtilities::SplitRectangle(
        std::vector<std::vector<CoordinatesArrayType>> &new_rectangles_master_parameter_space,
        std::vector<CoordinatesArrayType> rectangle_coordinates,
        double knot_line_position,
        bool is_vertical)
{
    // If the rectangle is not intersected by a knot line, return the complete rectangle
    if (!IsRectangleIntersectedByKnotLine(rectangle_coordinates, knot_line_position, is_vertical)){
        new_rectangles_master_parameter_space.push_back(rectangle_coordinates);
        return;
    }

    std::vector<CoordinatesArrayType> points_below_knot_line, points_above_knot_line;
    // std::vector<CoordinatesArrayType> intersections;
    // IndexType intersectionCount=0;

    // Categorize the coordinates of the given triangle
    categorizePoint(rectangle_coordinates[0], points_below_knot_line, points_above_knot_line, knot_line_position, is_vertical);
    categorizePoint(rectangle_coordinates[1], points_below_knot_line, points_above_knot_line, knot_line_position, is_vertical);
    categorizePoint(rectangle_coordinates[2], points_below_knot_line, points_above_knot_line, knot_line_position, is_vertical);
    categorizePoint(rectangle_coordinates[3], points_below_knot_line, points_above_knot_line, knot_line_position, is_vertical);

    if (is_vertical){
        CoordinatesArrayType intersection_point_1{knot_line_position, points_below_knot_line[0][1], 0.0};
        CoordinatesArrayType intersection_point_2{knot_line_position, points_below_knot_line[1][1], 0.0};

        // Append these intersection points to the vector of points below the knot line
        points_below_knot_line.push_back(intersection_point_1);
        points_below_knot_line.push_back(intersection_point_2);

        // Append these intersection points to the vector of points above the knot line
        points_above_knot_line.push_back(intersection_point_1);
        points_above_knot_line.push_back(intersection_point_2);
    } else {
        CoordinatesArrayType intersection_point_1{points_below_knot_line[0][0], knot_line_position, 0.0};
        CoordinatesArrayType intersection_point_2{points_below_knot_line[1][0], knot_line_position, 0.0};

        // Append these intersection points to the vector of points below the knot line
        points_below_knot_line.push_back(intersection_point_1);
        points_below_knot_line.push_back(intersection_point_2);

        // Append these intersection points to the vector of points above the knot line
        points_above_knot_line.push_back(intersection_point_1);
        points_above_knot_line.push_back(intersection_point_2);
    }

    std::vector<CoordinatesArrayType> new_rectangle_1 = OrderRectangleCoordinates(points_below_knot_line);
    std::vector<CoordinatesArrayType> new_rectangle_2 = OrderRectangleCoordinates(points_above_knot_line);

    // Append these newly created rectangles to the return object
    new_rectangles_master_parameter_space.push_back(new_rectangle_1);
    new_rectangles_master_parameter_space.push_back(new_rectangle_2);
}


bool IgaMappingIntersectionUtilities::IsRectangleIntersectedByKnotLine(
    const std::vector<CoordinatesArrayType>& rRectangleCoordinates,
    const double KnotLinePosition,
    const bool IsVertical)
{
    KRATOS_DEBUG_ERROR_IF(rRectangleCoordinates.size() != 4)
        << "Rectangle must have exactly 4 coordinates, but got "
        << rRectangleCoordinates.size() << std::endl;

    const double min_xi = std::min({
        rRectangleCoordinates[0][0],
        rRectangleCoordinates[1][0],
        rRectangleCoordinates[2][0],
        rRectangleCoordinates[3][0]
    });

    const double max_xi = std::max({
        rRectangleCoordinates[0][0],
        rRectangleCoordinates[1][0],
        rRectangleCoordinates[2][0],
        rRectangleCoordinates[3][0]
    });

    const double min_eta = std::min({
        rRectangleCoordinates[0][1],
        rRectangleCoordinates[1][1],
        rRectangleCoordinates[2][1],
        rRectangleCoordinates[3][1]
    });

    const double max_eta = std::max({
        rRectangleCoordinates[0][1],
        rRectangleCoordinates[1][1],
        rRectangleCoordinates[2][1],
        rRectangleCoordinates[3][1]
    });

    // Tolerance scaled with the size of the interval
    if (IsVertical) {
        const double length_xi = max_xi - min_xi;
        const double tolerance = 1e-12 * std::max(1.0, length_xi);

        // Degenerate rectangle in xi-direction
        if (length_xi <= tolerance) {
            return false;
        }

        return (KnotLinePosition > min_xi + tolerance) &&
               (KnotLinePosition < max_xi - tolerance);
    } else {
        const double length_eta = max_eta - min_eta;
        const double tolerance = 1e-12 * std::max(1.0, length_eta);

        // Degenerate rectangle in eta-direction
        if (length_eta <= tolerance) {
            return false;
        }

        return (KnotLinePosition > min_eta + tolerance) &&
               (KnotLinePosition < max_eta - tolerance);
    }
}

void IgaMappingIntersectionUtilities::categorizePoint(CoordinatesArrayType point,
                                                          std::vector<CoordinatesArrayType> &points_below,
                                                          std::vector<CoordinatesArrayType> &points_above,
                                                          double knot_line_position,
                                                          bool is_vertical)
{
    double tolerance = 1e-8;
    if ((is_vertical && point[0] < (knot_line_position - tolerance)) || (!is_vertical && point[1] < (knot_line_position - tolerance))){
        points_below.push_back(point);
    } else {
        if ((is_vertical && std::abs(point[0] - knot_line_position) < tolerance) || (!is_vertical && std::abs(point[1] - knot_line_position) < tolerance)){
            points_below.push_back(point);
        }
        points_above.push_back(point);
    }
}

std::vector<Kratos::IgaMappingIntersectionUtilities::CoordinatesArrayType> Kratos::IgaMappingIntersectionUtilities::OrderRectangleCoordinates(
    const std::vector<Kratos::IgaMappingIntersectionUtilities::CoordinatesArrayType>& rRectangleCoordinates)
{
    KRATOS_DEBUG_ERROR_IF(rRectangleCoordinates.size() != 4)
        << "Rectangle must have exactly 4 coordinates, but got "
        << rRectangleCoordinates.size() << std::endl;

    std::vector<CoordinatesArrayType> ordered_rectangle = rRectangleCoordinates;

    double center_x = 0.0;
    double center_y = 0.0;
    for (const auto& r_point : ordered_rectangle) {
        center_x += r_point[0];
        center_y += r_point[1];
    }
    center_x /= 4.0;
    center_y /= 4.0;

    std::sort(
        ordered_rectangle.begin(),
        ordered_rectangle.end(),
        [center_x, center_y](const CoordinatesArrayType& rA, const CoordinatesArrayType& rB) {
            const double angle_a = std::atan2(rA[1] - center_y, rA[0] - center_x);
            const double angle_b = std::atan2(rB[1] - center_y, rB[0] - center_x);
            return angle_a < angle_b;
        });

    return ordered_rectangle;
}

void IgaMappingIntersectionUtilities::WriteRectanglesToPythonPlot(
    const std::vector<std::vector<CoordinatesArrayType>>& rRectangles,
    const std::string& rFileName)
{
    std::ofstream out(rFileName);
    KRATOS_ERROR_IF_NOT(out.is_open())
        << "Could not open file: " << rFileName << std::endl;

    out << "import matplotlib.pyplot as plt\n";
    out << "rectangles = [\n";
    out << std::setprecision(16);

    for (IndexType i = 0; i < rRectangles.size(); ++i) {
        const auto& r_rect = rRectangles[i];

        KRATOS_ERROR_IF(r_rect.size() != 4)
            << "Rectangle " << i << " does not have 4 points, got "
            << r_rect.size() << std::endl;

        out << "    [";
        for (IndexType j = 0; j < 4; ++j) {
            out << "(" << r_rect[j][0] << ", " << r_rect[j][1] << ")";
            if (j != 3) {
                out << ", ";
            }
        }
        out << "]";
        if (i + 1 != rRectangles.size()) {
            out << ",";
        }
        out << "\n";
    }

    out << "]\n\n";

    out << "fig, ax = plt.subplots()\n";
    out << "for i, rect in enumerate(rectangles):\n";
    out << "    xs = [p[0] for p in rect] + [rect[0][0]]\n";
    out << "    ys = [p[1] for p in rect] + [rect[0][1]]\n";
    out << "    ax.plot(xs, ys, '-o')\n";
    out << "    cx = sum(p[0] for p in rect) / 4.0\n";
    out << "    cy = sum(p[1] for p in rect) / 4.0\n";
    out << "    ax.text(cx, cy, f'R{i}', fontsize=9)\n";
    out << "    for j, p in enumerate(rect):\n";
    out << "        ax.text(p[0], p[1], f'  {j}', fontsize=8)\n";
    out << "ax.set_aspect('equal')\n";
    out << "ax.grid(True)\n";
    out << "ax.set_xlabel('xi')\n";
    out << "ax.set_ylabel('eta')\n";
    out << "ax.set_title('Subdivided rectangles in parameter space')\n";
    out << "plt.show()\n";

    out.close();
}

} // namespace Kratos.
