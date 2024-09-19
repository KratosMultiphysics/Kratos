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
//                   Andrea Gorgi
//

// System includes

// External includes

// Project includes
#include "iga_mapping_intersection_utilities.h"

namespace Kratos
{

    void IgaMappingIntersectionUtilities::IgaCreateBrepCurveOnSurfaceCouplingGeometries(
        ModelPart &rModelPartDomainA,
        ModelPart &rModelPartDomainB,
        ModelPart &rModelPartResult,
        double Tolerance)
    {
        for (auto condition_a_itr = rModelPartDomainA.ConditionsBegin();
             condition_a_itr != rModelPartDomainA.ConditionsEnd();
             ++condition_a_itr)
        {
            // If domain B is weakly supported, the geometry is stored as a condition. On the other side, if it is strongly supported, the geometry is stored as a geometry
            if (rModelPartDomainB.NumberOfConditions() > 0){
                for (auto condition_b_itr = rModelPartDomainB.ConditionsBegin();
                    condition_b_itr != rModelPartDomainB.ConditionsEnd();
                    ++condition_b_itr)
                {
                    rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                        condition_a_itr->pGetGeometry(), condition_b_itr->pGetGeometry()));
                }
            } else{
                for (auto geometry_b_itr = rModelPartDomainB.GeometriesBegin();
                    geometry_b_itr != rModelPartDomainB.GeometriesEnd();
                    ++ geometry_b_itr)
                {
                    rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                        condition_a_itr->pGetGeometry(), rModelPartDomainB.pGetGeometry(geometry_b_itr->Id())));
                }
            }
        }
    }

    void IgaMappingIntersectionUtilities::CreateFEMIgaSurfaceCouplingGeometries(
        ModelPart &rModelPartDomainA,
        ModelPart &rModelPartDomainB,
        ModelPart &rModelPartResult,
        bool origin_is_iga)
    {
        // Get a pointer to the brep surface starting from the underlying element or condition (quadrature point)
        if (origin_is_iga == true)
        {
            GeometryType &brep_surface = rModelPartDomainA.ConditionsBegin()->pGetGeometry()->GetGeometryParent(0);

            std::vector<IndexType> patchs_id;

            // Create a vector with the ids of the different brep surfaces
            for (auto condition_itr = rModelPartDomainA.ConditionsBegin(); condition_itr != rModelPartDomainA.ConditionsEnd(); condition_itr++)
            {
                if (condition_itr == rModelPartDomainA.ConditionsBegin())
                {
                    patchs_id.push_back(condition_itr->pGetGeometry()->GetGeometryParent(0).Id());
                    continue;
                }

                std::ostringstream condition_name_stream;
                condition_itr->pGetGeometry()->PrintInfo(condition_name_stream);
                std::string condition_name = condition_name_stream.str();

                // We should avoid calling the GetGeometryParent(0) method for the Point load conditions!!
                if (condition_name == "a point in 3D space")
                {
                    continue;
                }

                IndexType current_patch_id = condition_itr->pGetGeometry()->GetGeometryParent(0).Id();

                auto find_patch_id_itr = std::find(patchs_id.begin(), patchs_id.end(), current_patch_id);

                if (find_patch_id_itr != patchs_id.end())
                {
                    continue; // The patch id is already inside the vector
                }
                else
                {
                    patchs_id.push_back(current_patch_id);
                }
            }

            KRATOS_WATCH(patchs_id)

            for (IndexType i = 0; i < patchs_id.size(); i++)
            {
                // Get a pointer to the brep surface geometry
                auto p_brep_surface = rModelPartDomainA.GetRootModelPart().pGetGeometry(patchs_id[i]);

                // Cast the nurbs surface pointer and brep surface to the derived class
                auto brep_surface_cast = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(p_brep_surface);

                // Creating one coupling geometry for each combination of finite element and patch surface
                /*for (auto element_b_itr = rModelPartDomainB.ElementsBegin();
                    element_b_itr != rModelPartDomainB.ElementsEnd();
                    ++ element_b_itr)
                {
                    rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                        p_brep_surface, element_b_itr->pGetGeometry()));
                }*/


                // Creating one coupling geometry for each combination of load condition and patch surface
                for (auto condition_b_itr = rModelPartDomainB.ConditionsBegin();
                     condition_b_itr != rModelPartDomainB.ConditionsEnd();
                     ++condition_b_itr)
                {
                    rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                        p_brep_surface, condition_b_itr->pGetGeometry()));
                }
            }
        }
        else
        {
            GeometryType &brep_surface = rModelPartDomainB.ConditionsBegin()->pGetGeometry()->GetGeometryParent(0);

            std::vector<IndexType> patchs_id;

            // Create a vector with the ids of the different brep surfaces
            for (auto condition_itr = rModelPartDomainB.ConditionsBegin(); condition_itr != rModelPartDomainB.ConditionsEnd(); condition_itr++)
            {
                if (condition_itr == rModelPartDomainB.ConditionsBegin())
                {
                    patchs_id.push_back(condition_itr->pGetGeometry()->GetGeometryParent(0).Id());
                    continue;
                }

                std::ostringstream condition_name_stream;
                condition_itr->pGetGeometry()->PrintInfo(condition_name_stream);
                std::string condition_name = condition_name_stream.str();

                // We should avoid calling the GetGeometryParent(0) method for the Point load conditions!!
                if (condition_name == "a point in 3D space")
                {
                    continue;
                }

                IndexType current_patch_id = condition_itr->pGetGeometry()->GetGeometryParent(0).Id();

                auto find_patch_id_itr = std::find(patchs_id.begin(), patchs_id.end(), current_patch_id);

                if (find_patch_id_itr != patchs_id.end())
                {
                    continue; // The patch id is already inside the vector
                }
                else
                {
                    patchs_id.push_back(current_patch_id);
                }
            }

            for (IndexType i = 0; i < patchs_id.size(); i++)
            {
                // Get a pointer to the brep surface geometry
                auto p_brep_surface = rModelPartDomainB.GetRootModelPart().pGetGeometry(patchs_id[i]);

                // Cast the nurbs surface pointer and brep surface to the derived class
                auto brep_surface_cast = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(p_brep_surface);

                // Creating one coupling geometry for each combination of finite element and patch surface
                /*for (auto element_a_itr = rModelPartDomainA.ElementsBegin();
                    element_a_itr != rModelPartDomainA.ElementsEnd();
                    ++ element_a_itr)
                {
                    rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                        p_brep_surface, element_a_itr->pGetGeometry()));
                }*/

                // Creating one coupling geometry for each combination of load condition and patch surface
                for (auto condition_b_itr = rModelPartDomainB.ConditionsBegin();
                     condition_b_itr != rModelPartDomainB.ConditionsEnd();
                     ++condition_b_itr)
                {
                    rModelPartResult.AddGeometry(Kratos::make_shared<CouplingGeometry<NodeType>>(
                        p_brep_surface, condition_b_itr->pGetGeometry()));
                }
            }
        }
    }

    void IgaMappingIntersectionUtilities::IgaCreateQuadraturePointsCoupling1DGeometries2D(
        ModelPart &rModelPartCoupling,
        double Tolerance)
    {
        const ModelPart &rParentModelPart = rModelPartCoupling.GetParentModelPart();

        for (auto geometry_itr = rModelPartCoupling.GeometriesBegin();
             geometry_itr != rModelPartCoupling.GeometriesEnd();
             ++geometry_itr)
        {
            std::ostringstream geometry_name_stream;
            geometry_itr->PrintInfo(geometry_name_stream);
            std::string geometry_name = geometry_name_stream.str();

            if (geometry_name == "Brep face curve"){continue;}

            IntegrationPointsArrayType integration_points;
            IntegrationInfo integration_info = geometry_itr->GetDefaultIntegrationInfo();

            geometry_itr->CreateIntegrationPoints(integration_points, integration_info);

            GeometriesArrayType master_and_slave_quadrature_points_geometries(integration_points.size()); // Vector of coupling geometries which stores the master and slave quedrature point geometries
            IndexType number_of_shape_functions_derivatives = 3;

            if (integration_points.size() != 0)
            {
                geometry_itr->CreateQuadraturePointGeometries(master_and_slave_quadrature_points_geometries, number_of_shape_functions_derivatives, integration_points, integration_info);

                const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
                                         ? 1
                                         : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;

                const SizeType IntegrationPointsPerSpan = integration_info.GetNumberOfIntegrationPointsPerSpan(0);

                for (IndexType i = 0; i < IntegrationPointsPerSpan; ++i)
                {
                    rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                        id + i, master_and_slave_quadrature_points_geometries(i)));
                }
            }
        }
    }

    /* Algorithm:
           - Check if the elements (knot spans) in the destination domain are trimmed or untrimmed
           - If they are trimmed, perform triangulation
           - Map to the physical space of the destination domain (rectangles if the domain is untrimmed or triangles if the domain is trimmed)
           - "Project" these geometrical figures to the origin domain
           - Clip the rectangles or triangles with the knot lines in the origin domain
           - Create quadrature point geometries in the origin domain
           - Map the integration points to the origin physical space
           - "Project" the integration points to the destination parameter space
           - Create quadrature point geometries in the destination domain
           - Create new "conditions" whose geometry is CouplingGeometry (origin and destination quadrature points)
   */

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
            if (geometry_itr->HasGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX))
            {
                continue;
            }

            // Extract the master and slave parts from the coupling geometry
            auto geom_master = geometry_itr->pGetGeometryPart(0); // IGA patch
            auto geom_slave = geometry_itr->pGetGeometryPart(1);  // IGA patch

            auto &r_geom_slave = geometry_itr->GetGeometryPart(1);

            auto geom_master_cast = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(geom_master);
            auto geom_slave_cast = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(geom_slave);

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
                    auto new_int_point = BrepTrimmingUtilities::ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);
                    if (!(int_point.x == new_int_point.x && int_point.y == new_int_point.y))
                    {
                        all_loops[0].push_back(new_int_point);
                        int_point.x = new_int_point.x;
                        int_point.y = new_int_point.y;
                    }
                }
            }

            // Loop over the "elements" in the destination domain
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

                    if (solution.size() > 0)
                    {
                        clip_area = std::abs(Clipper2Lib::Area(solution[0]));
                        for (IndexType k = 1; k < solution.size(); ++k)
                        {
                            clip_area -= std::abs(Clipper2Lib::Area(solution[k]));
                        }
                    }
                    else
                    {
                        continue; // If the solution size is equal to 0, continue to the next iteration
                    }

                    // Check if the element is untrimmed. If it is untrimmed, map the rectangle to the destination physical space and then project to the origin parameter space
                    if (std::abs(clip_area - span_area) < 1000)
                    {
                        // Create a vector with the coordinates of the rectangle in the destination physical space
                        std::vector<CoordinatesArrayType> rectangle_coordinates_xyz, rectangle_coordinates_master_parameter_space;

                        // Map the rectangle vertices from the parameter space to the physical space
                        for (IndexType k = i; k < i + 2; k++)
                        {
                            for (IndexType h = j; h < j + 2; h++)
                            {
                                CoordinatesArrayType vertex_coordinate_xyz, vertex_coordinate_slave_parameter_space{knot_vector_xi_slave_domain[k], knot_vector_eta_slave_domain[h], 0.0};

                                geom_slave->GlobalCoordinates(vertex_coordinate_xyz, vertex_coordinate_slave_parameter_space);

                                rectangle_coordinates_xyz.push_back(vertex_coordinate_xyz);
                            }
                        }

                        // Project the rectangle to the origin domain parameter space
                        for (IndexType k = 0; k < rectangle_coordinates_xyz.size(); k++)
                        {
                            CoordinatesArrayType vertex_coordinate_master_parameter_space;

                            geom_master->ProjectionPointGlobalToLocalSpace(rectangle_coordinates_xyz[k], vertex_coordinate_master_parameter_space, 1e-6);

                            rectangle_coordinates_master_parameter_space.push_back(vertex_coordinate_master_parameter_space);
                        }

                        // Check if the projected rectangles are cut by knot lines in the origin parameter space and subdivide them in new rectangles
                        std::vector<std::vector<CoordinatesArrayType>> obtained_rectangles;
                        IgaMappingIntersectionUtilities::SubdivideRectangleWithMasterKnotLines(rectangle_coordinates_master_parameter_space, geom_master, obtained_rectangles);

                        // Iterate over the newly created rectangles in the origin parameter space
                        for (IndexType k = 0; k < obtained_rectangles.size(); k++)
                        {
                            // Create integration points in the origin parameter space
                            IntegrationInfo master_integration_info = geom_master->GetDefaultIntegrationInfo();
                            IntegrationInfo slave_integration_info = geom_slave->GetDefaultIntegrationInfo();
                            const IndexType number_of_integration_points_master = master_integration_info.GetNumberOfIntegrationPointsPerSpan(0) * master_integration_info.GetNumberOfIntegrationPointsPerSpan(1);
                            const IndexType number_of_integration_points_slave = master_integration_info.GetNumberOfIntegrationPointsPerSpan(0) * master_integration_info.GetNumberOfIntegrationPointsPerSpan(1);
                            IndexType number_of_integration_points;

                            // Define the number of integration points to place in the origin domain as the biggest number between the origin and destination domains
                            if (number_of_integration_points_master >= number_of_integration_points_slave)
                            {
                                number_of_integration_points = number_of_integration_points_master;
                            }
                            else
                            {
                                number_of_integration_points = number_of_integration_points_slave;
                            }

                            // Create integration points and quadrature points containers for the master
                            IntegrationPointsArrayType integration_points_master(number_of_integration_points);
                            typename IntegrationPointsArrayType::iterator integration_point_master_iterator = integration_points_master.begin();
                            GeometriesArrayType quadrature_point_geometries_master(number_of_integration_points);

                            // Create integration and quadrature points containers for the slave
                            IntegrationPointsArrayType integration_points_slave(number_of_integration_points);
                            GeometriesArrayType quadrature_point_geometries_slave(number_of_integration_points);

                            if (obtained_rectangles[k][0][0] != obtained_rectangles[k][2][0])
                            {
                                // This is true when the physical space is not a rotation of the parameter space
                                IntegrationPointUtilities::IntegrationPoints2D(
                                    integration_point_master_iterator,
                                    master_integration_info.GetNumberOfIntegrationPointsPerSpan(0), master_integration_info.GetNumberOfIntegrationPointsPerSpan(1),
                                    obtained_rectangles[k][0][0], obtained_rectangles[k][2][0],
                                    obtained_rectangles[k][0][1], obtained_rectangles[k][1][1]);
                            }
                            else
                            {
                                IntegrationPointUtilities::IntegrationPoints2D(
                                    integration_point_master_iterator,
                                    master_integration_info.GetNumberOfIntegrationPointsPerSpan(0), master_integration_info.GetNumberOfIntegrationPointsPerSpan(1),
                                    obtained_rectangles[k][0][0], obtained_rectangles[k][1][0],
                                    obtained_rectangles[k][0][1], obtained_rectangles[k][2][1]);
                            }

                            // Create master quadrature points
                            geom_master->CreateQuadraturePointGeometries(quadrature_point_geometries_master, 3, integration_points_master, master_integration_info);

                            // Get the position of the master quadrature points in its physical space
                            CoordinatesArrayType master_quadrature_point_xyz = ZeroVector(3);
                            CoordinatesArrayType slave_quadrature_point_local_space = ZeroVector(3);
                            integration_points_slave = integration_points_master;

                            for (IndexType k = 0; k < quadrature_point_geometries_master.size(); k++)
                            {
                                // Get the position of the master i_th quadrature point in its physical space
                                master_quadrature_point_xyz = quadrature_point_geometries_master[k].Center();

                                // Obtain the local coordinate of the i_th master quadrature point in the slave local space
                                geom_slave->ProjectionPointGlobalToLocalSpace(master_quadrature_point_xyz, slave_quadrature_point_local_space, 1e-6);

                                // Modify the slave integration points container
                                integration_points_slave[k].X() = slave_quadrature_point_local_space[0];
                                integration_points_slave[k].Y() = slave_quadrature_point_local_space[1];
                            }

                            // Create slave quadrature points
                            geom_slave->CreateQuadraturePointGeometries(quadrature_point_geometries_slave, 3, integration_points_slave, slave_integration_info);

                            // Create new conditions in the coupling model part whose geometry is CouplingGeometry
                            const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
                                                     ? 1
                                                     : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;

                            for (IndexType i = 0; i < quadrature_point_geometries_master.size(); ++i)
                            {
                                rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                                    id + i, Kratos::make_shared<CouplingGeometry<Node>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i))));
                            }
                        }
                    }
                    else
                    {
                        std::vector<std::vector<CoordinatesArrayType>> triangles_destination_parameter_space; // Vector of triangles coming from the triangulation process (="triangles", but of different kind)
                        std::vector<Matrix> triangles;

                        // Do triangulation
                        BrepTrimmingUtilities::Triangulate_OPT(solution[0], triangles, factor);

                        for (IndexType u = 0; u < triangles.size(); ++u)
                        {
                            std::vector<CoordinatesArrayType> triangle; // One specific triangle coming from the triangulation process
                            // Extract the triangles from the Triangulation_OPT method
                            for (IndexType vertex_index = 0; vertex_index < 3; vertex_index++)
                            {
                                double coordinate_xi = triangles[u](vertex_index, 0);
                                double coordinate_eta = triangles[u](vertex_index, 1);
                                CoordinatesArrayType node_coordinate{coordinate_xi, coordinate_eta, 0.0};
                                triangle.push_back(node_coordinate);
                            }

                            // Order the triangles vertices counterclockwise
                            std::vector<CoordinatesArrayType> sorted_triangle;
                            sortVerticesCounterClockwise(triangle, sorted_triangle);

                            triangles_destination_parameter_space.push_back(sorted_triangle);
                        }

                        for (IndexType k = 0; k < triangles_destination_parameter_space.size(); k++)
                        {

                            // Create a vector with the coordinates of the triangle in the destination physical space
                            std::vector<CoordinatesArrayType> triangle_coordinates_xyz, triangle_coordinates_master_parameter_space;

                            // Map each triangle node from parameter to physical space
                            for (IndexType h = 0; h < 3; h++)
                            {
                                CoordinatesArrayType vertex_coordinate_xyz;
                                geom_slave->GlobalCoordinates(vertex_coordinate_xyz, triangles_destination_parameter_space[k][h]);
                                triangle_coordinates_xyz.push_back(vertex_coordinate_xyz);
                            }

                            // Project the triangle to the origin domain parameter space
                            for (IndexType k = 0; k < triangle_coordinates_xyz.size(); k++)
                            {
                                CoordinatesArrayType vertex_coordinate_master_parameter_space;

                                geom_master->ProjectionPointGlobalToLocalSpace(triangle_coordinates_xyz[k], vertex_coordinate_master_parameter_space, 1e-6);

                                triangle_coordinates_master_parameter_space.push_back(vertex_coordinate_master_parameter_space);
                            }

                            // Triangulation step: Intersect the triangles in the origin parameter space with the knot lines and subdivide them
                            std::vector<std::vector<CoordinatesArrayType>> obtained_triangles; // The outer vector contains the newly found triangles

                            // The triangulation process will create new triangles which are not intersected by the knot lines
                            IgaMappingIntersectionUtilities::Triangulation(triangle_coordinates_master_parameter_space, geom_master, obtained_triangles);

                            for (IndexType j = 0; j < obtained_triangles.size(); j++)
                            {
                                // Create integration points container for the master
                                IntegrationPointsArrayType integration_points_master(3);
                                typename IntegrationPointsArrayType::iterator integration_point_master_iterator = integration_points_master.begin();
                                GeometriesArrayType quadrature_point_geometries_master(3);

                                // Create integration points container for the slave
                                IntegrationPointsArrayType integration_points_slave(3);
                                GeometriesArrayType quadrature_point_geometries_slave(3);

                                // Coordinates of the triangle vertex in the master parameter space
                                double xi_0 = obtained_triangles[j][0][0];
                                double xi_1 = obtained_triangles[j][1][0];
                                double xi_2 = obtained_triangles[j][2][0];
                                double eta_0 = obtained_triangles[j][0][1];
                                double eta_1 = obtained_triangles[j][1][1];
                                double eta_2 = obtained_triangles[j][2][1];

                                // Creating the integration points
                                IntegrationPointUtilities::IntegrationPointsTriangle2D(integration_point_master_iterator, 1, xi_0, xi_1, xi_2, eta_0, eta_1, eta_2);

                                // Create master quadrature points
                                IntegrationInfo master_integration_info = geom_master->GetDefaultIntegrationInfo();
                                geom_master->CreateQuadraturePointGeometries(quadrature_point_geometries_master, 3, integration_points_master, master_integration_info);

                                // Get the position of the master quadrature points in its physical space
                                CoordinatesArrayType master_quadrature_point_xyz = ZeroVector(3);
                                CoordinatesArrayType slave_quadrature_point_local_space = ZeroVector(3);
                                integration_points_slave = integration_points_master;

                                for (IndexType k = 0; k < quadrature_point_geometries_master.size(); k++)
                                {
                                    // Get the position of the master i_th quadrature point in its physical space

                                    master_quadrature_point_xyz = quadrature_point_geometries_master[k].Center();

                                    // Obtain the local coordinate of the i_th master quadrature point in the slave local space
                                    geom_slave->ProjectionPointGlobalToLocalSpace(master_quadrature_point_xyz, slave_quadrature_point_local_space, 1e-6);

                                    // Modify the slave integration points container
                                    integration_points_slave[k].X() = slave_quadrature_point_local_space[0];
                                    integration_points_slave[k].Y() = slave_quadrature_point_local_space[1];
                                }

                                // Create slave quadrature points
                                IntegrationInfo slave_integration_info = geom_slave->GetDefaultIntegrationInfo();
                                geom_slave->CreateQuadraturePointGeometries(quadrature_point_geometries_slave, 3, integration_points_slave, slave_integration_info);

                                // add the quadrature point geometry conditions to the result model part
                                const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
                                                         ? 1
                                                         : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;

                                for (IndexType i = 0; i < quadrature_point_geometries_master.size(); ++i)
                                {
                                    rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                                        id + i, Kratos::make_shared<CouplingGeometry<Node>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i))));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void IgaMappingIntersectionUtilities::IgaFEMCreateQuadraturePointsCoupling2DGeometries3D(
        ModelPart &rModelPartCoupling,
        bool origin_is_iga)
    {
        const ModelPart &rParentModelPart = rModelPartCoupling.GetParentModelPart();

        // Iterate over the coupling geometries and create the quadrature points
        for (auto geometry_itr = rModelPartCoupling.GeometriesBegin();
             geometry_itr != rModelPartCoupling.GeometriesEnd();
             ++geometry_itr)
        {
            // Don´t enter the for loop if the geometry is a brep curve on surface or a brep surface (remember that in order to get the brep surface we have to include in the coupling model part the geometries of origin and destination)
            if (geometry_itr->HasGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX) )
            {
                continue;
            }

            auto geom_master = geometry_itr->pGetGeometryPart(0); // IGA Surface
            auto geom_slave = geometry_itr->pGetGeometryPart(1);  // Finite element

            KRATOS_WATCH(*geom_slave)
            KRATOS_WATCH(*geom_master)

            auto &r_geom_slave = geometry_itr->GetGeometryPart(1);

            auto geom_master_cast = dynamic_pointer_cast<BrepSurface<PointerVector<NodeType>, PointerVector<Point>>>(geom_master);
            auto geom_slave_cast = dynamic_pointer_cast<Triangle2D3<NodeType>>(geom_slave);

            GeometryPointerType master_nurbs_surface = geom_master_cast->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
            auto master_nurbs_surface_cast = dynamic_pointer_cast<NurbsSurfaceGeometry<3, PointerVector<Node>>>(master_nurbs_surface);
            
            CoordinatesArrayType local_parameter = ZeroVector(3);
            std::vector<std::vector<CoordinatesArrayType>> triangle_coordinates_parameter_space;
            std::vector<CoordinatesArrayType> triangle_coordinates;

            IndexType projection_is_successful_count = 0;
            std::vector<CoordinatesArrayType> projection_is_successful_nodes_coordinates_xyz;
            std::vector<CoordinatesArrayType> projection_is_not_successful_nodes_coordinates_xyz;

            // Create a vector which contains the points to do triangulation
            std::vector<CoordinatesArrayType> points_to_triangulate;

            for (IndexType i = 0; i < geom_slave->size(); i++)
            {
                CoordinatesArrayType node_coordinate_xyz = geom_slave->pGetPoint(i)->GetInitialPosition();

                KRATOS_WATCH(geom_slave->pGetPoint(i)->Id())

                IgaMappingIntersectionUtilities::FindInitialGuessNewtonRaphsonProjection(node_coordinate_xyz, geom_master, local_parameter);

                // Project each FEM node to the parameter space of the IGA Surface
                IndexType projection_is_successful = geom_master->ProjectionPointGlobalToLocalSpace(node_coordinate_xyz, local_parameter, 1e-6);

                if (projection_is_successful == 1)
                {
                    projection_is_successful_count += 1;
                    projection_is_successful_nodes_coordinates_xyz.push_back(geom_slave->pGetPoint(i)->GetInitialPosition());
                    points_to_triangulate.push_back(local_parameter);
                }
                else
                {
                    projection_is_not_successful_nodes_coordinates_xyz.push_back(geom_slave->pGetPoint(i)->GetInitialPosition());
                }
            }

            KRATOS_WATCH(geom_master->Id())
            KRATOS_WATCH(projection_is_successful_count)
            KRATOS_WATCH(points_to_triangulate)

            // If no node is projected inside the patch, continue to the next iteration. Otherwise, find the intersection between the triangle sides and the surface
            if (projection_is_successful_count == 0)
            {
                continue;
            }
            else if (projection_is_successful_count == 1)
            {   
                // If just one node is projected inside the surface, create one triangle

                // Check if the projection is on the boundary of the parameter space and, if so, continue to the next element
                if (AreProjectionsOnParameterSpaceBoundary(points_to_triangulate, *master_nurbs_surface_cast) == true) continue;

                std::vector<CoordinatesArrayType> triangle_sides_intersection_with_surface_patch_local_parameter;
                for (IndexType i = 0; i < 2; i++)
                {
                    CoordinatesArrayType intersection_point_local_space = ZeroVector(3);
                    IgaMappingIntersectionUtilities::FindIntersectionTriangleSideWithSurfacePatchBisectionMethod(geom_master, projection_is_successful_nodes_coordinates_xyz[0], projection_is_not_successful_nodes_coordinates_xyz[i], points_to_triangulate[0], intersection_point_local_space);
                    triangle_sides_intersection_with_surface_patch_local_parameter.push_back(intersection_point_local_space);
                    points_to_triangulate.push_back(intersection_point_local_space);
                }

                // If the intersection between the triangle and the patch surface is at just one point, ignore the element and go to the next triangle
                if (norm_2(triangle_sides_intersection_with_surface_patch_local_parameter[0] - triangle_sides_intersection_with_surface_patch_local_parameter[1]) < 1e-4)
                {
                    continue;
                }

                // Sort the vertices counter-clockwise before doing triangulation
                std::vector<CoordinatesArrayType> sorted_points_to_triangulate;
                sortVerticesCounterClockwise(points_to_triangulate, sorted_points_to_triangulate);

                // Form one triangle with the one point projected inside the patch and the two intersections
                triangle_coordinates_parameter_space.push_back(sorted_points_to_triangulate);
            }
            else if (projection_is_successful_count == 2)
            { 
                // If 2 nodes are projected inside the surface, create 2 triangles
                
                // Check if the projections are on the boundary of the parameter space and, if so, continue to the next element
                if (AreProjectionsOnParameterSpaceBoundary(points_to_triangulate, *master_nurbs_surface_cast) == true) continue;

                std::vector<CoordinatesArrayType> triangle_sides_intersection_with_surface_patch_local_parameter;
                for (IndexType i = 0; i < 2; i++)
                {
                    CoordinatesArrayType intersection_point_local_space = ZeroVector(3);
                    IgaMappingIntersectionUtilities::FindIntersectionTriangleSideWithSurfacePatchBisectionMethod(geom_master, projection_is_successful_nodes_coordinates_xyz[i], projection_is_not_successful_nodes_coordinates_xyz[0], points_to_triangulate[0],intersection_point_local_space);
                    if (norm_2(points_to_triangulate[i] - intersection_point_local_space) > 1e-4 && norm_2(points_to_triangulate.back() - intersection_point_local_space) > 1e-4)
                    {
                        points_to_triangulate.push_back(intersection_point_local_space);
                    }
                }

                KRATOS_WATCH(points_to_triangulate)
                // Sort the vertices counter-clockwise before doing triangulation
                std::vector<CoordinatesArrayType> sorted_points_to_triangulate;
                sortVerticesCounterClockwise(points_to_triangulate, sorted_points_to_triangulate);

                if (points_to_triangulate.size() == 3)
                {
                    // Create just one triangle
                    triangle_coordinates_parameter_space.push_back(sorted_points_to_triangulate);
                }
                else
                {
                    // Create two triangles
                    std::vector<CoordinatesArrayType> triangle_1, triangle_2;
                    triangle_1.insert(triangle_1.end(), {sorted_points_to_triangulate[0], sorted_points_to_triangulate[1], sorted_points_to_triangulate[2]});
                    triangle_2.insert(triangle_2.end(), {sorted_points_to_triangulate[0], sorted_points_to_triangulate[2], sorted_points_to_triangulate[3]});

                    // Insert the two triangles in one vector
                    triangle_coordinates_parameter_space.insert(triangle_coordinates_parameter_space.end(), {triangle_1, triangle_2});
                }
            }
            else if (projection_is_successful_count == 3)
            { // If the 3 nodes are projected inside the untrimmed surface, check if the triangles intersect the trimming curves and , if so, create new triangles

                // Get the trimming curves outer loop
                auto outer_loop_array = geom_master_cast->GetOuterLoops();

                Clipper2Lib::Paths64 all_loops(1), triangle(1), solution;
                const double factor = 1e-10;

                Clipper2Lib::Point64 int_point;
                int_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
                int_point.y = static_cast<cInt>(std::numeric_limits<int>::min());

                // Loop over the trimming curves and do tessellation
                for (IndexType j = 0; j < outer_loop_array[0].size(); ++j)
                {
                    CurveTessellation<PointerVector<Node>> curve_tesselation;
                    auto geometry_outer = *(outer_loop_array[0][j].get());
                    curve_tesselation.Tessellate(
                        geometry_outer, 0.001, 1, true);
                    auto tesselation = curve_tesselation.GetTessellation();
                    for (IndexType u = 0; u < tesselation.size(); ++u)
                    {
                        auto new_int_point = BrepTrimmingUtilities::ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);

                        // Add a certain tolerance to the curve
                        new_int_point.x += 100000.0;
                        new_int_point.y += 100000.0;

                        if (!(int_point.x == new_int_point.x && int_point.y == new_int_point.y))
                        {
                            all_loops[0].push_back(new_int_point);
                            int_point.x = new_int_point.x;
                            int_point.y = new_int_point.y;
                        }
                    }
                }

                // Build the triangle for clipper
                for (IndexType u = 0; u < 3; ++u)
                {
                    auto new_int_point = BrepTrimmingUtilities::ToIntPoint(points_to_triangulate[u][0], points_to_triangulate[u][1], factor);
                    if (!(int_point.x == new_int_point.x && int_point.y == new_int_point.y))
                    {
                        triangle[0].push_back(new_int_point);
                        int_point.x = new_int_point.x;
                        int_point.y = new_int_point.y;
                    }
                }

                KRATOS_WATCH(all_loops)

                // Intersect the triangle with the polygon resulting from the tessellation of the outer loop
                solution = Clipper2Lib::Intersect(all_loops, triangle, Clipper2Lib::FillRule::NonZero);

                // If an intersection is found, triangulate the clipped domain
                if (solution.size() != 0)
                {
                    std::vector<Matrix> triangles;
                    BrepTrimmingUtilities::Triangulate_OPT(solution[0], triangles, factor);

                    for (IndexType u = 0; u < triangles.size(); ++u)
                    {
                        // Extract the triangles from the Triangulation_OPT method
                        for (IndexType vertex_index = 0; vertex_index < 3; vertex_index++)
                        {
                            points_to_triangulate[vertex_index][0] = triangles[u](vertex_index, 0);
                            points_to_triangulate[vertex_index][1] = triangles[u](vertex_index, 1);
                            points_to_triangulate[vertex_index][2] = 0.0;
                        }

                        // Order the triangles vertices counterclockwise
                        std::vector<CoordinatesArrayType> sorted_points_to_triangulate;
                        sortVerticesCounterClockwise(points_to_triangulate, sorted_points_to_triangulate);

                        triangle_coordinates_parameter_space.push_back(sorted_points_to_triangulate);
                    }
                }
                else
                {
                   KRATOS_WATCH("I am here because the clipping is failing")
                   exit(0);
                   continue;
                }
            }

            KRATOS_WATCH(triangle_coordinates_parameter_space)

            for (IndexType i = 0; i < triangle_coordinates_parameter_space.size(); i++)
            {
                std::vector<std::vector<CoordinatesArrayType>> obtained_triangles; // The outer vector contains the newly found triangles

                // The triangulation process will create new triangles which are not intersected by the knot lines
                //IgaMappingIntersectionUtilities::Triangulation(triangle_coordinates_parameter_space[i], geom_master, obtained_triangles);
                obtained_triangles.push_back(triangle_coordinates_parameter_space[i]);

                for (IndexType j = 0; j < obtained_triangles.size(); j++)
                {
                    // Create integration points container for the master
                    IntegrationPointsArrayType integration_points_master(3);
                    typename IntegrationPointsArrayType::iterator integration_point_master_iterator = integration_points_master.begin();
                    GeometriesArrayType quadrature_point_geometries_master(3);

                    // Create integration points container for the slave
                    IntegrationPointsArrayType integration_points_slave(3);
                    GeometriesArrayType quadrature_point_geometries_slave(3);

                    // Coordinates of the triangle vertex in the master parameter space
                    double xi_0 = obtained_triangles[j][0][0];
                    double xi_1 = obtained_triangles[j][1][0];
                    double xi_2 = obtained_triangles[j][2][0];
                    double eta_0 = obtained_triangles[j][0][1];
                    double eta_1 = obtained_triangles[j][1][1];
                    double eta_2 = obtained_triangles[j][2][1];

                    // Creating the integration points
                    IntegrationPointUtilities::IntegrationPointsTriangle2D(integration_point_master_iterator, 1, xi_0, xi_1, xi_2, eta_0, eta_1, eta_2);

                    KRATOS_WATCH(integration_points_master)

                    // Create master quadrature points
                    IntegrationInfo master_integration_info = geom_master->GetDefaultIntegrationInfo();
                    geom_master->CreateQuadraturePointGeometries(quadrature_point_geometries_master, 3, integration_points_master, master_integration_info);

                    // Get the position of the master quadrature points in its physical space
                    CoordinatesArrayType master_quadrature_point_xyz = ZeroVector(3);
                    CoordinatesArrayType slave_quadrature_point_local_space = ZeroVector(3);
                    integration_points_slave = integration_points_master;

                    for (IndexType k = 0; k < quadrature_point_geometries_master.size(); k++)
                    {
                        // Get the position of the master i_th quadrature point in its physical space

                        master_quadrature_point_xyz = quadrature_point_geometries_master[k].Center();

                        // Obtain the local coordinate of the i_th master quadrature point in the slave local space
                        geom_slave->PointLocalCoordinates(slave_quadrature_point_local_space, master_quadrature_point_xyz);

                        // Modify the slave integration points container
                        integration_points_slave[k].X() = slave_quadrature_point_local_space[0];
                        integration_points_slave[k].Y() = slave_quadrature_point_local_space[1];
                    }

                    // Create slave quadrature points
                    CreateQuadraturePointsUtility<NodeType>::Create(
                        r_geom_slave, quadrature_point_geometries_slave, integration_points_slave, 1);

                    // add the quadrature point geometry conditions to the result model part
                    const IndexType id = (rParentModelPart.NumberOfConditions() == 0)
                                             ? 1
                                             : (rParentModelPart.ConditionsEnd() - 1)->Id() + 1;
                    if (origin_is_iga == true)
                    {
                        for (IndexType i = 0; i < quadrature_point_geometries_master.size(); ++i)
                        {
                            rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                                id + i, Kratos::make_shared<CouplingGeometry<Node>>(quadrature_point_geometries_master(i), quadrature_point_geometries_slave(i))));
                        }
                    }
                    else
                    {
                        rModelPartCoupling.AddCondition(Kratos::make_intrusive<Condition>(
                            id + i, Kratos::make_shared<CouplingGeometry<Node>>(quadrature_point_geometries_slave(i), quadrature_point_geometries_master(i))));
                    }
                }
            }
        }
        /*for (IndexType j = 0; j < nodes_to_add.size(); j++){
            rModelPartCoupling.GetSubModelPart("interface_destination").AddNode(nodes_to_add[j]);
        }*/
    }

    void IgaMappingIntersectionUtilities::FindInitialGuessNewtonRaphsonProjection(CoordinatesArrayType slave_element_node,
                                                          GeometryPointerType master_geometry,
                                                          CoordinatesArrayType &initial_guess)
    {
        std::vector<double> knot_vector_xi;
        std::vector<double> knot_vector_eta;

        // Points to evaluate
        std::vector<double> xi_coordinates;
        std::vector<double> eta_coordinates;

        // Obtain the knot vectors in the xi and eta directions
        master_geometry->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->SpansLocalSpace(knot_vector_xi, 0);
        master_geometry->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->SpansLocalSpace(knot_vector_eta, 1);

        // Construct the vector xi_coordinates and eta_coordinates 
        double xi_dist = (knot_vector_xi.back() - knot_vector_xi[0])/100;
        double eta_dist = (knot_vector_eta.back() - knot_vector_eta[0])/100;

        for(IndexType i = 0; i < 101; i++){
            xi_coordinates.push_back(knot_vector_xi[0]+(xi_dist*i));
        }

         for(IndexType j = 0; j < 101; j++){
            eta_coordinates.push_back(knot_vector_eta[0]+(eta_dist*j));
        }


        // Distance between slave_element_node and the position in the physial space of the considered point in the local space
        CoordinatesArrayType distance;

        CoordinatesArrayType local_parameter = ZeroVector(3);
        CoordinatesArrayType physical_space_position = ZeroVector(3);

        double min_distance = std::numeric_limits<double>::max();

        for(IndexType i = 0; i < xi_coordinates.size(); i++){
            for(IndexType j = 0; j < eta_coordinates.size(); j++){
                local_parameter[0] = xi_coordinates[i];
                local_parameter[1] = eta_coordinates[j];
                master_geometry->GlobalCoordinates(physical_space_position, local_parameter);
                
                 // Compute the distance vector
                CoordinatesArrayType distance_vector = physical_space_position - slave_element_node;

                // Compute the distance
                const double distance = norm_2(distance_vector);

                if (distance < min_distance){
                    min_distance = distance;
                    initial_guess = local_parameter;
                }
            }   
        }

        KRATOS_WATCH(min_distance)
        KRATOS_WATCH(initial_guess)
    }

    bool IgaMappingIntersectionUtilities::AreProjectionsOnParameterSpaceBoundary(std::vector<CoordinatesArrayType>& points_to_triangulate, NurbsSurfaceGeometry<3, PointerVector<Node>>& nurbs_surface)
    {
        NurbsInterval surface_interval_u = nurbs_surface.DomainIntervalU();
        NurbsInterval surface_interval_v = nurbs_surface.DomainIntervalU();

        // get the maximum and minimum values of these intervals
        double min_parameter_u = surface_interval_u.MinParameter();
        double max_parameter_u = surface_interval_u.MaxParameter();
        double min_parameter_v = surface_interval_v.MinParameter();
        double max_parameter_v = surface_interval_v.MaxParameter();

        if (points_to_triangulate.size() == 1){
            if (points_to_triangulate[0][0] == min_parameter_u || points_to_triangulate[0][0] == max_parameter_u) return true;
            if (points_to_triangulate[0][1] == min_parameter_v || points_to_triangulate[0][1] == max_parameter_v) return true;
        }
        else{
            if (points_to_triangulate[0][0] == min_parameter_u && points_to_triangulate[1][0] == min_parameter_u) return true;
            if (points_to_triangulate[0][0] == max_parameter_u && points_to_triangulate[1][0] == max_parameter_u) return true;
            if (points_to_triangulate[0][1] == min_parameter_v && points_to_triangulate[1][1] == min_parameter_v) return true;
            if (points_to_triangulate[0][1] == max_parameter_v && points_to_triangulate[1][1] == max_parameter_v) return true;
        }

        return false;
    }

    void IgaMappingIntersectionUtilities::Triangulation(
        std::vector<CoordinatesArrayType> original_triangle_coordinates,
        GeometryPointerType master_geometry,
        std::vector<std::vector<CoordinatesArrayType>> &new_triangles)
    {
        new_triangles.push_back(original_triangle_coordinates);

        std::vector<double> knot_vector_xi;
        std::vector<double> knot_vector_eta;

        // Obtain the knot vectors in the xi and eta directions
        master_geometry->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->SpansLocalSpace(knot_vector_xi, 0);
        master_geometry->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->SpansLocalSpace(knot_vector_eta, 1);

        // Check if both knot vectors size is 2 and if true, return the original triangle coordinates (no triangulation will de found)
        if (knot_vector_xi.size() == 2 && knot_vector_eta.size() == 2)
        {
            return;
        }

        std::vector<double> inner_knots_xi, inner_knots_eta;

        // Obtain the inner knots for the xi and eta direction
        if (knot_vector_xi.size() > 2)
        {
            inner_knots_xi = std::vector<double>(knot_vector_xi.begin() + 1, knot_vector_xi.end() - 1);
        }

        if (knot_vector_eta.size() > 2)
        {
            inner_knots_eta = std::vector<double>(knot_vector_eta.begin() + 1, knot_vector_eta.end() - 1);
        }

        for (IndexType k = 0; k < inner_knots_xi.size(); k++)
        {
            std::vector<std::vector<CoordinatesArrayType>> intermediate_triangles_xi;
            for (IndexType j = 0; j < new_triangles.size(); j++)
            {
                std::vector<std::vector<CoordinatesArrayType>> new_temp_triangles_xi;
                splitTriangle(new_temp_triangles_xi, new_triangles[j], inner_knots_xi[k], true);
                for (IndexType i = 0; i < new_temp_triangles_xi.size(); i++)
                {
                    intermediate_triangles_xi.push_back(new_temp_triangles_xi[i]);
                }
            }
            new_triangles = intermediate_triangles_xi;
        }

        for (IndexType k = 0; k < inner_knots_eta.size(); k++)
        {
            std::vector<std::vector<CoordinatesArrayType>> intermediate_triangles_eta;
            for (IndexType j = 0; j < new_triangles.size(); j++)
            {
                std::vector<std::vector<CoordinatesArrayType>> new_temp_triangles_eta;
                splitTriangle(new_temp_triangles_eta, new_triangles[j], inner_knots_eta[k], false);
                for (IndexType i = 0; i < new_temp_triangles_eta.size(); i++)
                {
                    intermediate_triangles_eta.push_back(new_temp_triangles_eta[i]);
                }
            }
            new_triangles = intermediate_triangles_eta;
        }
    }

    bool IgaMappingIntersectionUtilities::isTriangleIntersectedByKnotLine(
        std::vector<CoordinatesArrayType> triangle_coordinates,
        double knot_line_position,
        bool is_vertical)
    {
        double tolerance = 1e-8;
        if (is_vertical)
        {
            return ((triangle_coordinates[0][0] < (knot_line_position - tolerance) && triangle_coordinates[1][0] > (knot_line_position + tolerance)) ||
                    (triangle_coordinates[0][0] > (knot_line_position + tolerance) && triangle_coordinates[1][0] < (knot_line_position - tolerance)) ||
                    (triangle_coordinates[1][0] < (knot_line_position - tolerance) && triangle_coordinates[2][0] > (knot_line_position + tolerance)) ||
                    (triangle_coordinates[1][0] > (knot_line_position + tolerance) && triangle_coordinates[2][0] < (knot_line_position - tolerance)) ||
                    (triangle_coordinates[2][0] < (knot_line_position - tolerance) && triangle_coordinates[0][0] > (knot_line_position + tolerance)) ||
                    (triangle_coordinates[2][0] > (knot_line_position + tolerance) && triangle_coordinates[0][0] < (knot_line_position - tolerance)));
        }
        else
        {
            return ((triangle_coordinates[0][1] < (knot_line_position - tolerance) && triangle_coordinates[1][1] > (knot_line_position + tolerance)) ||
                    (triangle_coordinates[0][1] > (knot_line_position + tolerance) && triangle_coordinates[1][1] < (knot_line_position - tolerance)) ||
                    (triangle_coordinates[1][1] < (knot_line_position - tolerance) && triangle_coordinates[2][1] > (knot_line_position + tolerance)) ||
                    (triangle_coordinates[1][1] > (knot_line_position + tolerance) && triangle_coordinates[2][1] < (knot_line_position - tolerance)) ||
                    (triangle_coordinates[2][1] < (knot_line_position - tolerance) && triangle_coordinates[0][1] > (knot_line_position + tolerance)) ||
                    (triangle_coordinates[2][1] > (knot_line_position + tolerance) && triangle_coordinates[0][1] < (knot_line_position - tolerance)));
        }
    }

    void IgaMappingIntersectionUtilities::intersectionPoint(
        CoordinatesArrayType point_1,
        CoordinatesArrayType point_2,
        CoordinatesArrayType &intersection_point,
        double knot_line_position,
        bool is_vertical)
    {
        if (is_vertical)
        {
            double t = (knot_line_position - point_1[0]) / (point_2[0] - point_1[0]);
            intersection_point[0] = knot_line_position;
            intersection_point[1] = point_1[1] + t * (point_2[1] - point_1[1]);
            intersection_point[2] = 0.0;
        }
        else
        {
            double t = (knot_line_position - point_1[1]) / (point_2[1] - point_1[1]);
            intersection_point[0] = point_1[0] + t * (point_2[0] - point_1[0]);
            intersection_point[1] = knot_line_position;
            intersection_point[2] = 0.0;
        }
    }

    void IgaMappingIntersectionUtilities::categorizePoint(CoordinatesArrayType point,
                                                          std::vector<CoordinatesArrayType> &points_below,
                                                          std::vector<CoordinatesArrayType> &points_above,
                                                          double knot_line_position,
                                                          bool is_vertical)
    {
        double tolerance = 1e-8;
        if ((is_vertical && point[0] < (knot_line_position - tolerance)) || (!is_vertical && point[1] < (knot_line_position - tolerance)))
        {
            points_below.push_back(point);
        }
        else
        {
            if ((is_vertical && std::abs(point[0] - knot_line_position) < tolerance) || (!is_vertical && std::abs(point[1] - knot_line_position) < tolerance))
            {
                points_below.push_back(point);
            }
            points_above.push_back(point);
        }
    }

    void IgaMappingIntersectionUtilities::FindIntersectionTriangleSideWithSurfacePatchBisectionMethod(GeometryPointerType geom_master, CoordinatesArrayType point_inside, CoordinatesArrayType point_outside, CoordinatesArrayType initial_guess, CoordinatesArrayType &intersection_point)
    {
        CoordinatesArrayType mid_point_local_parameter = initial_guess;

        // Define a tolerance for the bisection method iterations
        double tolerance = 0;
        IndexType iterations_count = 0;

        CoordinatesArrayType new_mid_point;
        CoordinatesArrayType old_mid_point = ZeroVector(3);

        while (iterations_count < 100)
        {
            new_mid_point = 0.5 * (point_inside + point_outside);

            tolerance = norm_2(new_mid_point - old_mid_point);

            if (tolerance > 1e-8)
            {
                IndexType projection_is_successful = geom_master->ProjectionPointGlobalToLocalSpace(new_mid_point, mid_point_local_parameter, 1e-4);
                if (projection_is_successful == 0)
                {
                    point_outside = new_mid_point;
                }
                else
                {
                    point_inside = new_mid_point;
                }
                old_mid_point = new_mid_point;
                iterations_count++;
            }
            else
            {
                break;
            }
        }
        intersection_point = mid_point_local_parameter;
    }

    void IgaMappingIntersectionUtilities::sortVerticesCounterClockwise(std::vector<CoordinatesArrayType> &triangle_vertices, std::vector<CoordinatesArrayType> &sorted_triangle_vertices)
    {
        // Calculate the centroid of the triangle
        CoordinatesArrayType centroid = ZeroVector(3);

        for (IndexType i = 0; i < triangle_vertices.size(); i++)
        {
            centroid[0] += triangle_vertices[i][0];
            centroid[1] += triangle_vertices[i][1];
        }

        centroid[0] /= triangle_vertices.size();
        centroid[1] /= triangle_vertices.size();

        // Sort points counter-clockwise based on angle to centroid
        sorted_triangle_vertices = triangle_vertices;

        std::sort(sorted_triangle_vertices.begin(), sorted_triangle_vertices.end(), [&centroid](const CoordinatesArrayType &a, const CoordinatesArrayType &b)
                  { return atan2(a[1] - centroid[1], a[0] - centroid[0]) < atan2(b[1] - centroid[1], b[0] - centroid[0]); });
    }

    void IgaMappingIntersectionUtilities::splitTriangle(
        std::vector<std::vector<CoordinatesArrayType>> &new_triangles,
        std::vector<CoordinatesArrayType> triangle_coordinates,
        double knot_line_position,
        bool is_vertical)
    {
        // If the triangle is not intersected by a knot line, return the complete triangle
        if (!isTriangleIntersectedByKnotLine(triangle_coordinates, knot_line_position, is_vertical))
        {
            new_triangles.push_back(triangle_coordinates);
            return;
        }

        std::vector<CoordinatesArrayType> points_below, points_above;
        std::vector<CoordinatesArrayType> intersections;
        IndexType intersectionCount = 0;

        // Categorize the coordinates of the given triangle
        categorizePoint(triangle_coordinates[0], points_below, points_above, knot_line_position, is_vertical);
        categorizePoint(triangle_coordinates[1], points_below, points_above, knot_line_position, is_vertical);
        categorizePoint(triangle_coordinates[2], points_below, points_above, knot_line_position, is_vertical);

        double tolerance = 1e-8;

        if (is_vertical)
        {
            if ((triangle_coordinates[0][0] < (knot_line_position - tolerance) && triangle_coordinates[1][0] > (knot_line_position + tolerance)) ||
                (triangle_coordinates[0][0] > (knot_line_position + tolerance) && triangle_coordinates[1][0] < (knot_line_position - tolerance)))
            {
                CoordinatesArrayType intersection_point;
                intersectionPoint(triangle_coordinates[0], triangle_coordinates[1], intersection_point, knot_line_position, is_vertical);
                intersections.push_back(intersection_point);
                intersectionCount += 1;
            }
            if ((triangle_coordinates[1][0] < (knot_line_position - tolerance) && triangle_coordinates[2][0] > (knot_line_position + tolerance)) ||
                (triangle_coordinates[1][0] > (knot_line_position + tolerance) && triangle_coordinates[2][0] < (knot_line_position - tolerance)))
            {
                CoordinatesArrayType intersection_point;
                intersectionPoint(triangle_coordinates[1], triangle_coordinates[2], intersection_point, knot_line_position, is_vertical);
                intersections.push_back(intersection_point);
                intersectionCount += 1;
            }
            if ((triangle_coordinates[2][0] < (knot_line_position - tolerance) && triangle_coordinates[0][0] > (knot_line_position + tolerance)) ||
                (triangle_coordinates[2][0] > (knot_line_position + tolerance) && triangle_coordinates[0][0] < (knot_line_position - tolerance)))
            {
                CoordinatesArrayType intersection_point;
                intersectionPoint(triangle_coordinates[2], triangle_coordinates[0], intersection_point, knot_line_position, is_vertical);
                intersections.push_back(intersection_point);
                intersectionCount += 1;
            }
        }
        else
        {
            if ((triangle_coordinates[0][1] < (knot_line_position - tolerance) && triangle_coordinates[1][1] > (knot_line_position + tolerance)) ||
                (triangle_coordinates[0][1] > (knot_line_position + tolerance) && triangle_coordinates[1][1] < (knot_line_position - tolerance)))
            {
                CoordinatesArrayType intersection_point;
                intersectionPoint(triangle_coordinates[0], triangle_coordinates[1], intersection_point, knot_line_position, is_vertical);
                intersections.push_back(intersection_point);
                intersectionCount += 1;
            }
            if ((triangle_coordinates[1][1] < (knot_line_position - tolerance) && triangle_coordinates[2][1] > (knot_line_position + tolerance)) ||
                (triangle_coordinates[1][1] > (knot_line_position + tolerance) && triangle_coordinates[2][1] < (knot_line_position - tolerance)))
            {
                CoordinatesArrayType intersection_point;
                intersectionPoint(triangle_coordinates[1], triangle_coordinates[2], intersection_point, knot_line_position, is_vertical);
                intersections.push_back(intersection_point);
                intersectionCount += 1;
            }
            if ((triangle_coordinates[2][1] < (knot_line_position - tolerance) && triangle_coordinates[0][1] > (knot_line_position + tolerance)) ||
                (triangle_coordinates[2][1] > (knot_line_position + tolerance) && triangle_coordinates[0][1] < (knot_line_position - tolerance)))
            {
                CoordinatesArrayType intersection_point;
                intersectionPoint(triangle_coordinates[2], triangle_coordinates[0], intersection_point, knot_line_position, is_vertical);
                intersections.push_back(intersection_point);
                intersectionCount += 1;
            }
        }

        if (intersectionCount == 1)
        {
            points_below.push_back(intersections[0]);
            points_above.push_back(intersections[0]);
        }
        else if (intersectionCount == 2)
        {
            points_below.push_back(intersections[0]);
            points_above.push_back(intersections[0]);
            points_below.push_back(intersections[1]);
            points_above.push_back(intersections[1]);
        }

        if (points_below.size() == 3)
        { // Create just one new triangle
            std::vector<CoordinatesArrayType> new_triangle_below;
            std::vector<CoordinatesArrayType> sorted_points_below;
            sortVerticesCounterClockwise(points_below, sorted_points_below);

            new_triangle_below.push_back(sorted_points_below[0]);
            new_triangle_below.push_back(sorted_points_below[1]);
            new_triangle_below.push_back(sorted_points_below[2]);
            new_triangles.push_back(new_triangle_below);
        }
        else if (points_below.size() > 3)
        { // Create two new triangles
            std::vector<CoordinatesArrayType> new_triangle_below_1, new_triangle_below_2;
            std::vector<CoordinatesArrayType> sorted_points_below;
            sortVerticesCounterClockwise(points_below, sorted_points_below);

            // Create the first triangle
            new_triangle_below_1.push_back(sorted_points_below[0]);
            new_triangle_below_1.push_back(sorted_points_below[1]);
            new_triangle_below_1.push_back(sorted_points_below[2]);
            new_triangles.push_back(new_triangle_below_1);

            // Create the second triangle
            new_triangle_below_2.push_back(sorted_points_below[2]);
            new_triangle_below_2.push_back(sorted_points_below[3]);
            new_triangle_below_2.push_back(sorted_points_below[0]);
            new_triangles.push_back(new_triangle_below_2);
        }

        if (points_above.size() == 3)
        {
            std::vector<CoordinatesArrayType> new_triangle_above;
            std::vector<CoordinatesArrayType> sorted_points_above;
            sortVerticesCounterClockwise(points_above, sorted_points_above);

            new_triangle_above.push_back(sorted_points_above[0]);
            new_triangle_above.push_back(sorted_points_above[1]);
            new_triangle_above.push_back(sorted_points_above[2]);
            new_triangles.push_back(new_triangle_above);
        }
        else if (points_above.size() > 3)
        { // Create two new triangles
            std::vector<CoordinatesArrayType> new_triangle_above_1, new_triangle_above_2;
            std::vector<CoordinatesArrayType> sorted_points_above;
            sortVerticesCounterClockwise(points_above, sorted_points_above);

            // Create the first triangle
            new_triangle_above_1.push_back(sorted_points_above[0]);
            new_triangle_above_1.push_back(sorted_points_above[1]);
            new_triangle_above_1.push_back(sorted_points_above[2]);
            new_triangles.push_back(new_triangle_above_1);

            // Create the second triangle
            new_triangle_above_2.push_back(sorted_points_above[2]);
            new_triangle_above_2.push_back(sorted_points_above[3]);
            new_triangle_above_2.push_back(sorted_points_above[0]);
            new_triangles.push_back(new_triangle_above_2);
        }
    }

    void IgaMappingIntersectionUtilities::SubdivideRectangleWithMasterKnotLines(
        std::vector<CoordinatesArrayType> original_rectangle_coordinates,
        GeometryPointerType master_geometry,
        std::vector<std::vector<CoordinatesArrayType>> &new_rectangles)
    {
        new_rectangles.push_back(original_rectangle_coordinates);

        std::vector<double> knot_vector_xi;
        std::vector<double> knot_vector_eta;

        // Obtain the knot vectors in the xi and eta directions
        master_geometry->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->SpansLocalSpace(knot_vector_xi, 0);
        master_geometry->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX)->SpansLocalSpace(knot_vector_eta, 1);

        // Check if both knot vectors size is 2 and if true, return the original rectangle coordinates (no new rectangles will de found)
        if (knot_vector_xi.size() == 2 && knot_vector_eta.size() == 2)
        {
            return;
        }

        std::vector<double> inner_knots_xi, inner_knots_eta;

        // Obtain the inner knots for the xi and eta direction
        if (knot_vector_xi.size() > 2)
        {
            inner_knots_xi = std::vector<double>(knot_vector_xi.begin() + 1, knot_vector_xi.end() - 1);
        }

        if (knot_vector_eta.size() > 2)
        {
            inner_knots_eta = std::vector<double>(knot_vector_eta.begin() + 1, knot_vector_eta.end() - 1);
        }

        for (IndexType k = 0; k < inner_knots_xi.size(); k++)
        {
            std::vector<std::vector<CoordinatesArrayType>> intermediate_rectangles_xi;
            for (IndexType j = 0; j < new_rectangles.size(); j++)
            {
                std::vector<std::vector<CoordinatesArrayType>> new_temp_rectangles_xi;
                SplitRectangle(new_temp_rectangles_xi, new_rectangles[j], inner_knots_xi[k], true);
                for (IndexType i = 0; i < new_temp_rectangles_xi.size(); i++)
                {
                    intermediate_rectangles_xi.push_back(new_temp_rectangles_xi[i]);
                }
            }
            new_rectangles = intermediate_rectangles_xi;
        }

        for (IndexType k = 0; k < inner_knots_eta.size(); k++)
        {
            std::vector<std::vector<CoordinatesArrayType>> intermediate_rectangles_eta;
            for (IndexType j = 0; j < new_rectangles.size(); j++)
            {
                std::vector<std::vector<CoordinatesArrayType>> new_temp_rectangles_eta;
                SplitRectangle(new_temp_rectangles_eta, new_rectangles[j], inner_knots_eta[k], false);
                for (IndexType i = 0; i < new_temp_rectangles_eta.size(); i++)
                {
                    intermediate_rectangles_eta.push_back(new_temp_rectangles_eta[i]);
                }
            }
            new_rectangles = intermediate_rectangles_eta;
        }
    }

    void IgaMappingIntersectionUtilities::SplitRectangle(
        std::vector<std::vector<CoordinatesArrayType>> &new_rectangles,
        std::vector<CoordinatesArrayType> rectangle_coordinates,
        double knot_line_position,
        bool is_vertical)
    {
        // If the rectangle is not intersected by a knot line, return the complete rectangle
        if (!IsRectangleIntersectedByKnotLine(rectangle_coordinates, knot_line_position, is_vertical))
        {
            new_rectangles.push_back(rectangle_coordinates);
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

        if (is_vertical)
        {
            CoordinatesArrayType intersection_point_1{knot_line_position, points_below_knot_line[0][1], 0.0};
            CoordinatesArrayType intersection_point_2{knot_line_position, points_below_knot_line[1][1], 0.0};

            // Append these intersection points to the vector of points below the knot line
            points_below_knot_line.push_back(intersection_point_1);
            points_below_knot_line.push_back(intersection_point_2);

            // Append these intersection points to the vector of points above the knot line
            points_above_knot_line.push_back(intersection_point_1);
            points_above_knot_line.push_back(intersection_point_2);
        }
        else
        {
            CoordinatesArrayType intersection_point_1{points_below_knot_line[0][0], knot_line_position, 0.0};
            CoordinatesArrayType intersection_point_2{points_below_knot_line[1][0], knot_line_position, 0.0};

            // Append these intersection points to the vector of points below the knot line
            points_below_knot_line.push_back(intersection_point_1);
            points_below_knot_line.push_back(intersection_point_2);

            // Append these intersection points to the vector of points above the knot line
            points_above_knot_line.push_back(intersection_point_1);
            points_above_knot_line.push_back(intersection_point_2);
        }

        std::vector<CoordinatesArrayType> new_rectangle_1, new_rectangle_2;

        // Create new rectangles
        new_rectangle_1 = points_below_knot_line;
        new_rectangle_2 = points_above_knot_line;

        // Append these newly created rectangles to the return object
        new_rectangles.push_back(new_rectangle_1);
        new_rectangles.push_back(new_rectangle_2);
    }

    bool IgaMappingIntersectionUtilities::IsRectangleIntersectedByKnotLine(
        std::vector<CoordinatesArrayType> rectangle_coordinates,
        double knot_line_position,
        bool is_vertical)
    {
        // Get the bounding box of the rectangle
        double minXi = std::min({rectangle_coordinates[0][0], rectangle_coordinates[1][0], rectangle_coordinates[2][0], rectangle_coordinates[3][0]});
        double maxXi = std::max({rectangle_coordinates[0][0], rectangle_coordinates[1][0], rectangle_coordinates[2][0], rectangle_coordinates[3][0]});
        double minEta = std::min({rectangle_coordinates[0][1], rectangle_coordinates[1][1], rectangle_coordinates[2][1], rectangle_coordinates[3][1]});
        double maxEta = std::max({rectangle_coordinates[0][1], rectangle_coordinates[1][1], rectangle_coordinates[2][1], rectangle_coordinates[3][1]});

        if (is_vertical)
        {
            if (knot_line_position >= minXi && knot_line_position <= maxXi)
            {
                return true;
            }
        }
        else
        {
            if (knot_line_position >= minEta && knot_line_position <= maxEta)
            {
                return true;
            }
        }
    }

} // namespace Kratos.
