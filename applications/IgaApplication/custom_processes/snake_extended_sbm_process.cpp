//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "snake_extended_sbm_process.h"
#include "iga_application_variables.h"

namespace Kratos
{

SnakeExtendedSbmProcess::SnakeExtendedSbmProcess(
    Model& rModel, Parameters ThisParameters) : 
    SnakeSbmProcess(rModel, ThisParameters)
{

}

void SnakeExtendedSbmProcess::CreateSbmExtendedGeometries()
{
    FindClosestTruePointToSurrogateVertex();
}

void SnakeExtendedSbmProcess::FindClosestTruePointToSurrogateVertex()
{
    // get the data
    mEchoLevel = mThisParameters["echo_level"].GetInt();
    std::string iga_model_part_name = mThisParameters["model_part_name"].GetString();
    std::string skin_model_part_name = mThisParameters["skin_model_part_name"].GetString();
    std::string skin_model_part_inner_initial_name = mThisParameters["skin_model_part_inner_initial_name"].GetString();
    
    // Loop over the surrogate vertex and call the search in radius or something like that
    // TODO: Outer boyndary
    const bool is_inner = true;
    
    std::string surrogate_sub_model_part_name; 
    std::string skin_sub_model_part_name; 
    surrogate_sub_model_part_name = "surrogate_inner";
    skin_sub_model_part_name = "inner";

    //FIXME:
    ModelPart& r_extended_sbm_sub_model_part = mpIgaModelPart->CreateSubModelPart("extended_sbm");
    ModelPart& r_interface_extended_sbm_sub_model_part = mpIgaModelPart->CreateSubModelPart("interface_extended_sbm");

    ModelPart& r_skin_sub_model_part = mpSkinModelPart->GetSubModelPart(skin_sub_model_part_name);
    ModelPart& r_surrogate_sub_model_part = mpIgaModelPart->GetSubModelPart(surrogate_sub_model_part_name);

    // Create the testBins for the search in radius
    PointVector points;
    for (auto &i_cond : r_skin_sub_model_part.Conditions()) {
        points.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry().Center().X(), i_cond.GetGeometry().Center().Y(), i_cond.GetGeometry().Center().Z())));
    }
    // Get the mesh sizes from the surrogate model part
    const Vector& knot_span_sizes = r_surrogate_sub_model_part.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);

    double knot_span_reference_size = knot_span_sizes[0];
    if (knot_span_sizes[1] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[1];}
    if (knot_span_sizes.size() > 2) {if (knot_span_sizes[2] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[2];}}

    const int domain_size = mpIgaModelPart->GetProcessInfo()[DOMAIN_SIZE];
    double search_radius;
    if (domain_size == 2) {
        search_radius = 2*std::sqrt(2.0) * knot_span_reference_size;
    } else {
        KRATOS_ERROR << "This method is only implemented for 2D (DOMAIN_SIZE == 2). "
                    << "Current DOMAIN_SIZE: " << domain_size << std::endl;
    }

    DynamicBins testBins(points.begin(), points.end());

    // Maximum number of results to be found in the search in radius
    const int number_of_results = 1e6; 

    ModelPart::NodesContainerType::ContainerType results(number_of_results);
    std::vector<double> list_of_distances(number_of_results);

    auto p_surface = mpIgaModelPart->pGetGeometry(1);
    IndexType id_brep_curve_on_surface = (mpIgaModelPart->GeometriesEnd()-1)->Id() + 1;
    
    typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceType;
    auto p_nurbs_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(
                            p_surface->pGetGeometryPart(Geometry<typename PointerVector<NodeType>::value_type>::BACKGROUND_GEOMETRY_INDEX));
    IntegrationInfo surface_integration_info = p_nurbs_surface->GetDefaultIntegrationInfo();

    // std::string condition_name = "ExtendedSbmSolidCondition";
    std::string condition_name = "ExtendedSbmLoadSolidCondition";
    std::string interface_condition_name = "ExtendedSbmSolidInterfaceCondition";
    std::string element_name = "ExtendedSbmSolidElement";
    // std::string element_name = "SolidElement";

    // Loop over the nodes of the surrogate sub model part
    IndexType iel = 1;
    SizeType number_of_shape_functions_derivatives = 7;

    IndexType first_condition_id = r_surrogate_sub_model_part.pGetElement(iel)->GetGeometry()[0].Id();
    IndexType last_condition_id = r_surrogate_sub_model_part.pGetElement(iel)->GetGeometry()[1].Id();

    SizeType size_surrogate_loop = last_condition_id - first_condition_id + 1;

    for (SizeType j = 0; j < size_surrogate_loop; ++j) {
        auto p_brep_geometry = mpIgaModelPart->pGetGeometry(6 + j);
        auto p_brep_curve_on_surface_surrogate1_surrogate2 = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep_geometry);

        KRATOS_ERROR_IF(!p_brep_curve_on_surface_surrogate1_surrogate2) <<  ":::[SnakeExtendedSbmProcess]::: the geometry with id " << p_brep_curve_on_surface_surrogate1_surrogate2->Id() 
                                            << " is not a BrepCurveOnSurfaceType." << std::endl;

        NurbsInterval brep_domain_interval = p_brep_curve_on_surface_surrogate1_surrogate2->DomainInterval();
        CoordinatesArrayType surrogate_vertex_1 = ZeroVector(3); 
        CoordinatesArrayType surrogate_vertex_1_local_coords = ZeroVector(3);
        CoordinatesArrayType surrogate_vertex_2 = ZeroVector(3); 
        CoordinatesArrayType surrogate_vertex_2_local_coords = ZeroVector(3);
        surrogate_vertex_1_local_coords[0] = brep_domain_interval.GetT0();
        surrogate_vertex_2_local_coords[0] = brep_domain_interval.GetT1();

        p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_1, surrogate_vertex_1_local_coords);
        p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_2, surrogate_vertex_2_local_coords);

        // retrieve middle point of the brep
        CoordinatesArrayType surrogate_middle_point = ZeroVector(3); 
        CoordinatesArrayType surrogate_middle_point_local_coords = ZeroVector(3);
        surrogate_middle_point_local_coords[0] = 0.5 * (surrogate_vertex_1_local_coords[0] + surrogate_vertex_2_local_coords[0]);
        p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_middle_point, surrogate_middle_point_local_coords);

        IntegrationPoint<1> integration_point(surrogate_middle_point_local_coords[0]);
        IntegrationPointsArrayType surrogate_integration_points_list;
        surrogate_integration_points_list.push_back(integration_point);

        IntegrationInfo integration_info = p_brep_curve_on_surface_surrogate1_surrogate2->GetDefaultIntegrationInfo();
        GeometriesArrayType quadrature_point_list;
        p_brep_curve_on_surface_surrogate1_surrogate2->CreateQuadraturePointGeometries(quadrature_point_list, number_of_shape_functions_derivatives, 
                                                            surrogate_integration_points_list, integration_info);

        GeometryType::Pointer surrogate_brep_middle_geometry = quadrature_point_list(0);


        // STORE THE SURROGATE BREP MIDDLE GEOMETRY FOR THE LATERAL BREPS
        const double tol = 1.0e-12;
        bool check_cond_1 = false;
        bool check_cond_2 = false;
        for (auto& r_surrogate_node : r_surrogate_sub_model_part.Nodes())
        {
            // Check if this node coincides with surrogate_vertex_1
            if (norm_2(r_surrogate_node.Coordinates() - surrogate_vertex_1) < tol)
            {
                // Direct reference to the stored neighbour geometries
                auto& r_neighbour_geometries =
                    r_surrogate_node.GetValue(NEIGHBOUR_GEOMETRIES);

                // Append the new middle geometry; no SetValue call required
                r_neighbour_geometries.push_back(surrogate_brep_middle_geometry);

                check_cond_1 = true; // Only one node should satisfy the condition
            }
            else if (norm_2(r_surrogate_node.Coordinates() - surrogate_vertex_2) < tol)
            {
                // Direct reference to the stored neighbour geometries
                auto& r_neighbour_geometries =
                    r_surrogate_node.GetValue(NEIGHBOUR_GEOMETRIES);

                // Append the new middle geometry; no SetValue call required
                r_neighbour_geometries.push_back(surrogate_brep_middle_geometry);

                check_cond_2 = true; // Only one node should satisfy the condition
            }
            if (check_cond_1 && check_cond_2) {
                // If both conditions are satisfied, break the loop
                break;
            }
        }

        // search the projection of the first vertex
        DynamicBinsPointerType p_surrogate_vertex_1 = DynamicBinsPointerType(new PointType(1, surrogate_vertex_1[0], surrogate_vertex_1[1], surrogate_vertex_1[2]));
        // Apply the SearchInRadius
        SizeType obtained_results = testBins.SearchInRadius(*p_surrogate_vertex_1, search_radius, results.begin(), list_of_distances.begin(), number_of_results);
        double minimum_distance = 1e14;
        // Find the nearest node
        IndexType nearest_node_id;
        for (IndexType k = 0; k < obtained_results; k++) {
            double current_distance = list_of_distances[k];   
            if (current_distance < minimum_distance) { 
                minimum_distance = current_distance;
                nearest_node_id = k;
            }
        }
        KRATOS_ERROR_IF(obtained_results == 0) << "::[SnakeSbmProcess]:: Zero points found in serch for projection of point: " <<
            p_surrogate_vertex_1 << std::endl;
        
        const IndexType id_closest_true_node = results[nearest_node_id]->Id();
        const auto skin_vertex_1 = r_skin_sub_model_part.pGetNode(id_closest_true_node);
        
        // search the projection of the second vertex
        DynamicBinsPointerType p_surrogate_vertex_2 = DynamicBinsPointerType(new PointType(1, surrogate_vertex_2[0], surrogate_vertex_2[1], surrogate_vertex_2[2]));
        
        // Apply the SearchInRadius
        SizeType obtained_results_2 = testBins.SearchInRadius(*p_surrogate_vertex_2, search_radius, results.begin(), list_of_distances.begin(), number_of_results);
        
        minimum_distance = 1e14;
        // Find the nearest node
        nearest_node_id;
        for (IndexType k = 0; k < obtained_results_2; k++) {
            double current_distance = list_of_distances[k];   
            if (current_distance < minimum_distance) { 
                minimum_distance = current_distance;
                nearest_node_id = k;
            }
        }
        KRATOS_ERROR_IF(obtained_results_2 == 0) << "::[SnakeSbmProcess]:: Zero points found in serch for projection of point: " <<
            p_surrogate_vertex_2 << std::endl;
        
        const IndexType id_closest_true_node_2 = results[nearest_node_id]->Id();
        const auto skin_vertex_2 = r_skin_sub_model_part.pGetNode(id_closest_true_node_2);
    
        Vector active_range_knot_vector = ZeroVector(2);
        active_range_knot_vector[0] = 0;
        active_range_knot_vector[1] = 1;
        NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);

        // surrogate_1 - skin_1
        // create the brep connecting vertex and closest true point
        Point surrogate_1(surrogate_vertex_1);
        Point skin_1(*skin_vertex_1);
        
        Node::Pointer p_surrogate1_brep_point = Node::Pointer(new Node(1, surrogate_1));
        Node::Pointer p_skin1_brep_point = Node::Pointer(new Node(2, skin_1));

        auto p_nurbs_curve_surrogate1_skin1 = this->CreateBrepCurve(p_surrogate1_brep_point, p_skin1_brep_point, active_range_knot_vector);
        auto p_brep_curve_surrogate1_skin1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_skin1);      
                                                
        // surrogate_2 - skin_2
        // create the brep connecting vertex and closest true point
        Point surrogate_2(surrogate_vertex_2);
        Point skin_2(*skin_vertex_2);
        
        Node::Pointer p_surrogate2_brep_point = Node::Pointer(new Node(1, surrogate_2));
        Node::Pointer p_skin2_brep_point = Node::Pointer(new Node(2, skin_2));

        auto p_nurbs_curve_surrogate2_skin2 = this->CreateBrepCurve(p_surrogate2_brep_point, p_skin2_brep_point, active_range_knot_vector);
        auto p_brep_curve_surrogate2_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate2_skin2);      

        // skin_1 - skin_2
        // create the brep connecting vertex and closest true point

        IndexType id_skin_middle_point = 0.5*(id_closest_true_node + id_closest_true_node_2); 

        if (id_closest_true_node < id_closest_true_node_2) 
                id_skin_middle_point = 0.5*(r_skin_sub_model_part.NumberOfNodes() + id_closest_true_node + id_closest_true_node_2);

                if (id_skin_middle_point > r_skin_sub_model_part.Nodes().back().Id()) 
                    id_skin_middle_point -= r_skin_sub_model_part.NumberOfNodes();
                else if (id_skin_middle_point > r_skin_sub_model_part.Nodes().front().Id())
                    id_skin_middle_point += r_skin_sub_model_part.NumberOfNodes();


        if (id_skin_middle_point > r_skin_sub_model_part.Nodes().back().Id()) {
            id_skin_middle_point = id_skin_middle_point - r_skin_sub_model_part.Nodes().back().Id();
        }
        Node::Pointer p_skin_middle_brep_point = mpSkinModelPart->pGetNode(id_skin_middle_point);

        PointerVector<Node> ctrl_pts;
        Vector               knots;
        Vector               weights;

        BuildParabolicNurbsData(p_skin1_brep_point, p_skin_middle_brep_point, p_skin2_brep_point,
                                ctrl_pts, knots, weights, 0.5);

        const double polynomial_degree = 2;
        typename NurbsCurveGeometry<3, PointerVector<Node>>::Pointer p_nurbs_curve_skin1_skin2(
                new NurbsCurveGeometry<3, PointerVector<Node>>(
                    ctrl_pts,
                    polynomial_degree,
                    knots,
                    weights)); 

        // auto p_nurbs_curve_skin1_skin2 = this->CreateBrepCurve(p_skin1_brep_point, p_skin2_brep_point, active_range_knot_vector);
        auto p_brep_curve_skin1_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_skin1_skin2);      
    
        IntegrationPointsArrayType brep_integration_points_list_skin1_skin2;
        GeometriesArrayType brep_quadrature_point_list_skin1_skin2;

        p_brep_curve_skin1_skin2->CreateIntegrationPoints(brep_integration_points_list_skin1_skin2, integration_info);

        const double p_brep_curve_skin1_skin2_length = (p_brep_curve_skin1_skin2->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX))->Length();
        for (auto& integration_point : brep_integration_points_list_skin1_skin2) {
            integration_point.SetWeight(integration_point.Weight() * p_brep_curve_skin1_skin2_length);
        }

        p_brep_curve_skin1_skin2->CreateQuadraturePointGeometries(brep_quadrature_point_list_skin1_skin2, number_of_shape_functions_derivatives, 
                                                                  brep_integration_points_list_skin1_skin2, integration_info);
        
        SizeType id = 1;
        if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
            id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;
        
        std::vector<Geometry<Node>::Pointer> neighbour_geometries_skin1_skin2;
        neighbour_geometries_skin1_skin2.push_back(surrogate_brep_middle_geometry);
        
        this->CreateConditions(
            brep_quadrature_point_list_skin1_skin2.ptr_begin(), brep_quadrature_point_list_skin1_skin2.ptr_end(),
            r_extended_sbm_sub_model_part, condition_name, id, PropertiesPointerType(), knot_span_sizes, neighbour_geometries_skin1_skin2);

        
        // surrogate_1 - surrogate_2
        // create ONLY the brep connecting the two vertices
        auto p_nurbs_curve_surrogate1_surrogate2 = this->CreateBrepCurve(p_surrogate1_brep_point, p_surrogate2_brep_point, active_range_knot_vector);
        auto p_brep_curve_surrogate1_surrogate2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_surrogate2);      


        // Creation integration points on the cut elements

        IntegrationPointsArrayType surface_integration_points = CreateCoonsPatchGaussPoints(
                            3, /*Order*/
                            *p_brep_curve_surrogate1_surrogate2,   // B0
                            *p_brep_curve_surrogate1_skin1,       // L0
                            *p_brep_curve_surrogate2_skin2,       // L1
                            *p_brep_curve_skin1_skin2,            // B1
                            surrogate_1,  // P00
                            skin_1,       // P01
                            surrogate_2,  // P10
                            skin_2);      // P11
        
        GeometriesArrayType surface_quadrature_point_list;

        p_nurbs_surface->CreateQuadraturePointGeometries(surface_quadrature_point_list, number_of_shape_functions_derivatives, 
                                                        surface_integration_points, surface_integration_info);

        IndexType id_element = 1;
        if (mpIgaModelPart->GetRootModelPart().Elements().size() > 0)
            id_element = mpIgaModelPart->GetRootModelPart().Elements().back().Id() + 1;

        this->CreateElements(
            surface_quadrature_point_list.ptr_begin(), surface_quadrature_point_list.ptr_end(),
            *mpIgaModelPart, element_name, id_element, PropertiesPointerType(), neighbour_geometries_skin1_skin2);
    }
    
    //---------------------------------------------------------------------------
    bool is_entering = false;
    for (IndexType i_cond_id = first_condition_id; i_cond_id <= last_condition_id; ++i_cond_id)
    {
        is_entering = !is_entering;

        const auto& surrogate_condition = r_surrogate_sub_model_part.pGetCondition(i_cond_id);

        const auto& p_surrogate_1 = surrogate_condition->GetGeometry()(0);
        const auto& p_surrogate_2 = surrogate_condition->GetGeometry()(1); 
        
        const auto surrogate_middle_point = surrogate_condition->GetGeometry().Center();

        const bool is_first_surrogate_already_computed = (p_surrogate_1->GetValue(ACTIVATION_LEVEL) == 1); 
        const bool is_second_surrogate_already_computed = (p_surrogate_2->GetValue(ACTIVATION_LEVEL) == 1);

        // Project the surrogate vertices to the skin boundary
        // search the projection of the first vertex
        if (!is_first_surrogate_already_computed) {
            p_surrogate_1->SetValue(ACTIVATION_LEVEL, 1.0); 
            DynamicBinsPointerType p_surrogate_vertex_1 = DynamicBinsPointerType(new PointType(1, p_surrogate_1->X(), p_surrogate_1->Y(), p_surrogate_1->Z()));
            // Apply the SearchInRadius
            SizeType obtained_results = testBins.SearchInRadius(*p_surrogate_vertex_1, search_radius, results.begin(), list_of_distances.begin(), number_of_results);
            double minimum_distance = 1e14;
            // Find the nearest node
            IndexType nearest_node_id;
            for (IndexType k = 0; k < obtained_results; k++) {
                double current_distance = list_of_distances[k];   
                if (current_distance < minimum_distance) { 
                    minimum_distance = current_distance;
                    nearest_node_id = k;
                }
            }
            KRATOS_ERROR_IF(obtained_results == 0) << "::[SnakeSbmProcess]:: Zero points found in serch for projection of point: " <<
                p_surrogate_vertex_1 << std::endl;
            
            const IndexType id_closest_true_node = results[nearest_node_id]->Id();
            const auto skin_vertex_1 = r_skin_sub_model_part.pGetNode(id_closest_true_node);


            Vector active_range_knot_vector = ZeroVector(2);
            active_range_knot_vector[0] = 0;
            active_range_knot_vector[1] = 1;
            NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);

            // surrogate_1 - skin_1
            // create the brep connecting vertex and closest true point
            Point surrogate_1(*p_surrogate_1);
            Point skin_1(*skin_vertex_1);
            
            Node::Pointer p_first_point = Node::Pointer(new Node(1, surrogate_1));
            Node::Pointer p_second_point = Node::Pointer(new Node(2, skin_1));

            if (is_entering) {
                // change the order to preserve the anticlockwise orientation
                Node::Pointer p_temp_pointer = p_first_point;
                p_first_point = p_second_point;
                p_second_point = p_temp_pointer;
            } 

            auto p_nurbs_curve_surrogate1_skin1 = this->CreateBrepCurve(p_first_point, p_second_point, active_range_knot_vector);
            auto p_brep_curve_surrogate1_skin1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_skin1);      
        
            IntegrationInfo brep_integration_info_surrogate1_skin1 = p_brep_curve_surrogate1_skin1->GetDefaultIntegrationInfo();

            //FIXME: 
            const int p = 3;
            brep_integration_info_surrogate1_skin1.SetNumberOfIntegrationPointsPerSpan(0,2*p+1);

            IntegrationPointsArrayType brep_integration_points_list_surrogate1_skin1;
            GeometriesArrayType brep_quadrature_point_list_surrogate1_skin1;

            p_brep_curve_surrogate1_skin1->CreateIntegrationPoints(brep_integration_points_list_surrogate1_skin1, brep_integration_info_surrogate1_skin1);

            const double brep_curve_surrogate1_skin1_length = norm_2(skin_1 - surrogate_1);
            for (auto& integration_point : brep_integration_points_list_surrogate1_skin1) {
                integration_point.SetWeight(integration_point.Weight() * brep_curve_surrogate1_skin1_length);
            }

            p_brep_curve_surrogate1_skin1->CreateQuadraturePointGeometries(brep_quadrature_point_list_surrogate1_skin1, number_of_shape_functions_derivatives, 
                                                                brep_integration_points_list_surrogate1_skin1, brep_integration_info_surrogate1_skin1);
            
            SizeType id = 1;
            if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
                id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;

            auto& neighbour_geometries = p_surrogate_1->GetValue(NEIGHBOUR_GEOMETRIES);

            if (norm_2(surrogate_middle_point-neighbour_geometries[0]->Center()) < 1e-12) {
                // neighbour_geometries already disposed in the correct way
            } else if (norm_2(surrogate_middle_point-neighbour_geometries[1]->Center()) < 1e-12) {
                // neighbour_geometries not disposed in the correct way: reverse the entries
                const auto& temp_geom = neighbour_geometries[0];
                neighbour_geometries[0] = neighbour_geometries[1];
                neighbour_geometries[1] = temp_geom;

                p_surrogate_1->SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries);
            } else {
                KRATOS_ERROR << "::[SnakeSbmProcess]:: The surrogate middle point is not close to any of the neighbour geometries." << std::endl;
            }
            this->CreateConditions(
                brep_quadrature_point_list_surrogate1_skin1.ptr_begin(), brep_quadrature_point_list_surrogate1_skin1.ptr_end(),
                r_interface_extended_sbm_sub_model_part, interface_condition_name, id, PropertiesPointerType(), knot_span_sizes, neighbour_geometries);
        }

        if (!is_second_surrogate_already_computed) {
            p_surrogate_2->SetValue(ACTIVATION_LEVEL, 1.0); 
            DynamicBinsPointerType p_surrogate_vertex_2 = DynamicBinsPointerType(new PointType(1, p_surrogate_2->X(), p_surrogate_2->Y(), p_surrogate_2->Z()));
            // Apply the SearchInRadius
            SizeType obtained_results = testBins.SearchInRadius(*p_surrogate_vertex_2, search_radius, results.begin(), list_of_distances.begin(), number_of_results);
            double minimum_distance = 1e14;
            // Find the nearest node
            IndexType nearest_node_id;
            for (IndexType k = 0; k < obtained_results; k++) {
                double current_distance = list_of_distances[k];   
                if (current_distance < minimum_distance) { 
                    minimum_distance = current_distance;
                    nearest_node_id = k;
                }
            }
            KRATOS_ERROR_IF(obtained_results == 0) << "::[SnakeSbmProcess]:: Zero points found in serch for projection of point: " <<
                p_surrogate_vertex_2 << std::endl;
            
            const IndexType id_closest_true_node = results[nearest_node_id]->Id();
            const auto skin_vertex_2 = r_skin_sub_model_part.pGetNode(id_closest_true_node);


            Vector active_range_knot_vector = ZeroVector(2);
            active_range_knot_vector[0] = 0;
            active_range_knot_vector[1] = 1;
            NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);

            // surrogate_2 - skin_2
            // create the brep connecting vertex and closest true point
            Point surrogate_2(*p_surrogate_2);
            Point skin_2(*skin_vertex_2);

            Node::Pointer p_first_point = Node::Pointer(new Node(1, surrogate_2));
            Node::Pointer p_second_point = Node::Pointer(new Node(2, skin_2));

            if (!is_entering) {
                // change the order to preserve the anticlockwise orientation
                Node::Pointer p_temp_pointer = p_first_point;
                p_first_point = p_second_point;
                p_second_point = p_temp_pointer;
            } 

            // change the order to preserve the anticlockwise orientation
            auto p_nurbs_curve_surrogate2_skin2 = this->CreateBrepCurve(p_first_point, p_second_point, active_range_knot_vector);
            auto p_brep_curve_surrogate2_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate2_skin2);      
        
            IntegrationInfo brep_integration_info_surrogate2_skin2 = p_brep_curve_surrogate2_skin2->GetDefaultIntegrationInfo();

            //FIXME: 
            const int p = 3;
            brep_integration_info_surrogate2_skin2.SetNumberOfIntegrationPointsPerSpan(0,2*p+1);

            IntegrationPointsArrayType brep_integration_points_list_surrogate2_skin2;
            GeometriesArrayType brep_quadrature_point_list_surrogate2_skin2;

            p_brep_curve_surrogate2_skin2->CreateIntegrationPoints(brep_integration_points_list_surrogate2_skin2, brep_integration_info_surrogate2_skin2);

            const double brep_curve_surrogate2_skin2_length = norm_2(skin_2 - surrogate_2);
            for (auto& integration_point : brep_integration_points_list_surrogate2_skin2) {
                integration_point.SetWeight(integration_point.Weight() * brep_curve_surrogate2_skin2_length);
            }

            p_brep_curve_surrogate2_skin2->CreateQuadraturePointGeometries(brep_quadrature_point_list_surrogate2_skin2, number_of_shape_functions_derivatives, 
                                                                brep_integration_points_list_surrogate2_skin2, brep_integration_info_surrogate2_skin2);

                      
            SizeType id = 1;
            if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
                id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;

            
            auto& neighbour_geometries = p_surrogate_2->GetValue(NEIGHBOUR_GEOMETRIES);

            if (norm_2(surrogate_middle_point-neighbour_geometries[0]->Center()) < 1e-12) {
                // neighbour_geometries already disposed in the correct way
                
            } else if (norm_2(surrogate_middle_point-neighbour_geometries[1]->Center()) < 1e-12) {
                // neighbour_geometries not disposed in the correct way: reverse the entries
                const auto& temp_geom = neighbour_geometries[0];
                neighbour_geometries[0] = neighbour_geometries[1];
                neighbour_geometries[1] = temp_geom;

                p_surrogate_2->SetValue(NEIGHBOUR_GEOMETRIES, neighbour_geometries);
            } else {
                KRATOS_ERROR << "::[SnakeSbmProcess]:: The surrogate middle point is not close to any of the neighbour geometries." << std::endl;
            }

            this->CreateConditions(
                brep_quadrature_point_list_surrogate2_skin2.ptr_begin(), brep_quadrature_point_list_surrogate2_skin2.ptr_end(),
                r_interface_extended_sbm_sub_model_part, interface_condition_name, id, PropertiesPointerType(), knot_span_sizes, neighbour_geometries);

        }
    }
        
    
}


bool SnakeExtendedSbmProcess::ProjectToSkinBoundary(
        const ModelPart* pSkinModelPart,
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rProjectedPoint,
        CoordinatesArrayType& rProjectedPointLocal,
        std::vector<array_1d<double, 3>>& rCurveDerivatives,
        int nInitialGuesses)
{
    rProjectedPoint = ZeroVector(3);
    rCurveDerivatives.resize(3);
    bool is_projected_at_least_once = false;
    double best_distance = 1e12;
    std::vector<array_1d<double, 3>> best_curve_derivatives(2, ZeroVector(3));
    std::string best_layer_name = "";

    for (auto &i_curve : pSkinModelPart->Geometries())
    {   
        int nurbs_curve_id = i_curve.Id();
        auto p_nurbs_curve_geometry = pSkinModelPart->pGetGeometry(nurbs_curve_id);
        auto nurbs_curve_geometry = std::dynamic_pointer_cast<NurbsCurveGeometryType>(p_nurbs_curve_geometry);
        KRATOS_ERROR_IF(!nurbs_curve_geometry) <<  ":::[IgaContactProcessSbm]::: the geometry with id " << nurbs_curve_id 
                            << " is not a NurbsCurveGeometryType." << std::endl;

        const double t0 = nurbs_curve_geometry->DomainInterval().GetT0();
        const double t1 = nurbs_curve_geometry->DomainInterval().GetT1();

        for (int i_guess = 0; i_guess < nInitialGuesses; ++i_guess) {
            CoordinatesArrayType projected_point_local = ZeroVector(3);
            CoordinatesArrayType projected_point = ZeroVector(3);
            std::vector<array_1d<double, 3>> curve_derivatives(2, ZeroVector(3));

            projected_point_local[0] = t0 + (t1 - t0) * double(i_guess) / (nInitialGuesses - 1);

            bool is_projected = nurbs_curve_geometry->ProjectionPointGlobalToLocalSpace(rPoint, projected_point_local, 1e-13);

            if (!is_projected) continue;

            nurbs_curve_geometry->GlobalCoordinates(projected_point, projected_point_local);

            double curr_distance = norm_2(rPoint - projected_point);

            if (curr_distance < best_distance) {
                best_distance = curr_distance;
                rProjectedPoint = projected_point;
                nurbs_curve_geometry->GlobalSpaceDerivatives(best_curve_derivatives, projected_point_local, 2);
                rProjectedPointLocal = projected_point_local;
                is_projected_at_least_once = true;
                best_layer_name = i_curve.GetValue(IDENTIFIER);


            }
        }
    }

    rCurveDerivatives = best_curve_derivatives;

    if (!is_projected_at_least_once)
    {
        KRATOS_WARNING("::[IgaContactProcessSbm]:: no projection found on the skin boundary")  
                    << " for the point: " << rPoint << std::endl;
    }

    return is_projected_at_least_once;
}


void SnakeExtendedSbmProcess::CreateConditions(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    std::string& rConditionName,
    SizeType& rIdCounter,
    PropertiesPointerType pProperties,
    const Vector KnotSpanSizes,
    const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries) const
{
    const Condition& reference_condition = KratosComponents<Condition>::Get(rConditionName);

    ModelPart::ConditionsContainerType new_condition_list;

    KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
        << "Creating conditions of type " << rConditionName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    IndexType geometry_count = 0;
    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
    {
        new_condition_list.push_back(
            reference_condition.Create(rIdCounter, (*it), pProperties));
        
        // Set knot span sizes to the condition
        new_condition_list.GetContainer()[geometry_count]->SetValue(KNOT_SPAN_SIZES, KnotSpanSizes);

        new_condition_list.GetContainer()[geometry_count]->SetValue(NEIGHBOUR_GEOMETRIES, pSurrogateReferenceGeometries);

        rIdCounter++;
        geometry_count++;
    }

    rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
}


void SnakeExtendedSbmProcess::CreateElements(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    std::string& rElementName,
    SizeType& rIdCounter,
    PropertiesPointerType pProperties,
    const std::vector<Geometry<Node>::Pointer> &pSurrogateReferenceGeometries) const
{
    KRATOS_ERROR_IF(!KratosComponents<Element>::Has(rElementName))
        << rElementName << " not registered." << std::endl;

    const Element& rReferenceElement = KratosComponents<Element>::Get(rElementName);

    ElementsContainerType new_element_list;

    KRATOS_INFO_IF("CreateElements", mEchoLevel > 2)
        << "Creating elements of type " << rElementName
        << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

    SizeType num_elements = std::distance(rGeometriesBegin, rGeometriesEnd);
    new_element_list.reserve(num_elements);
    IndexType geometry_count = 0;

    for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
    {
        new_element_list.push_back(
            rReferenceElement.Create(rIdCounter, (*it), pProperties));

        new_element_list.GetContainer()[geometry_count]->SetValue(NEIGHBOUR_GEOMETRIES, pSurrogateReferenceGeometries);
        rIdCounter++;
        geometry_count++;
    }

    rModelPart.AddElements(new_element_list.begin(), new_element_list.end());
}


// ------------------------------------------------------------------
// 1-D Gauss–Legendre on [0,1]
// ------------------------------------------------------------------
void SnakeExtendedSbmProcess::GaussLegendreOnUnitInterval(
    const std::size_t      Order,
    std::vector<double>&   rXi,
    std::vector<double>&   rWeight)
{
    KRATOS_ERROR_IF(Order < 1 || Order > 5)
        << "Gauss order " << Order << " not implemented (1…5 supported)." << std::endl;

    //  nodi/mesi in [-1,1]  ------------------------------
    static const double sXi  [5][5] = {
        {  0.0,                                   0.0, 0.0, 0.0, 0.0 }, // n=1
        { -0.5773502691896257,  +0.5773502691896257,                    }, // n=2
        {  0.0, -0.7745966692414834, +0.7745966692414834,               }, // n=3
        { -0.3399810435848563, +0.3399810435848563,
          -0.8611363115940526, +0.8611363115940526,                     }, // n=4
        {  0.0,
          -0.5384693101056831, +0.5384693101056831,
          -0.9061798459386640, +0.9061798459386640 }                      // n=5
    };

    static const double sW   [5][5] = {
        {  2.0,                                   0.0, 0.0, 0.0, 0.0 }, // n=1
        {  1.0,  1.0,                                                }, // n=2
        {  0.8888888888888888, 0.5555555555555556, 0.5555555555555556 }, // n=3
        {  0.6521451548625461, 0.6521451548625461,
           0.3478548451374539, 0.3478548451374539                     }, // n=4
        {  0.5688888888888889,
           0.4786286704993665, 0.4786286704993665,
           0.2369268850561891, 0.2369268850561891 }                     // n=5
    };

    const std::size_t n = Order;

    rXi.resize(n);
    rWeight.resize(n);

    for (std::size_t k=0; k<n; ++k)
    {
        // mappa nodo da [-1,1] → [0,1]
        rXi[k]     = 0.5 * (sXi[n-1][k] + 1.0);
        // riscalo peso: w' = w / 2
        rWeight[k] = 0.5 *  sW [n-1][k];
    }
}

// ------------------------------------------------------------------
// Global point on Brep curve
// ------------------------------------------------------------------
array_1d<double,3> SnakeExtendedSbmProcess::GlobalPoint(
    const GeometryType& rCurve,
    const double         T)
{
    CoordinatesArrayType local(3,0.0);
    local[0] = T;

    array_1d<double,3> P;
    rCurve.GlobalCoordinates(P, local);
    return P;
}

// ------------------------------------------------------------------
// Coons patch mapping X(ξ,η)
// ------------------------------------------------------------------
array_1d<double,3> SnakeExtendedSbmProcess::CoonsPoint(
    const double                  Xi,
    const double                  Eta,
    const GeometryType&          rB0,
    const GeometryType&          rL0,
    const GeometryType&          rL1,
    const GeometryType&          rB1,
    const array_1d<double,3>&     rP00,
    const array_1d<double,3>&     rP01,
    const array_1d<double,3>&     rP10,
    const array_1d<double,3>&     rP11)
{
    const array_1d<double,3> B0 = GlobalPoint(rB0, Xi);
    const array_1d<double,3> B1 = GlobalPoint(rB1, Xi);
    const array_1d<double,3> L0 = GlobalPoint(rL0, Eta);
    const array_1d<double,3> L1 = GlobalPoint(rL1, Eta);

    const double om_xi  = 1.0 - Xi;
    const double om_eta = 1.0 - Eta;

    return  om_xi * L0 + Xi * L1
          + om_eta * B0 + Eta * B1
          - om_xi * om_eta * rP00
          - Xi    * om_eta * rP10
          - Xi    * Eta    * rP11
          - om_xi * Eta    * rP01;
}

// ------------------------------------------------------------------
// Finite-difference derivative ∂X/∂ξ or ∂X/∂η
// ------------------------------------------------------------------
array_1d<double,3> SnakeExtendedSbmProcess::CoonsDerivativeFD(
    const double                  Xi,
    const double                  Eta,
    const bool                    WithRespectToXi,
    const GeometryType&          rB0,
    const GeometryType&          rL0,
    const GeometryType&          rL1,
    const GeometryType&          rB1,
    const array_1d<double,3>&     rP00,
    const array_1d<double,3>&     rP01,
    const array_1d<double,3>&     rP10,
    const array_1d<double,3>&     rP11,
    const double                  Step)
{
    if (WithRespectToXi) {
        return 0.5/Step *
               ( CoonsPoint(Xi+Step, Eta, rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11) -
                 CoonsPoint(Xi-Step, Eta, rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11) );
    } else {
        return 0.5/Step *
               ( CoonsPoint(Xi, Eta+Step, rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11) -
                 CoonsPoint(Xi, Eta-Step, rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11) );
    }
}

// ------------------------------------------------------------------
// Build Gauss points on the curved quadrilateral (tensor-product)
// ------------------------------------------------------------------
SnakeExtendedSbmProcess::IntegrationPointsArrayType
SnakeExtendedSbmProcess::CreateCoonsPatchGaussPoints(
    const std::size_t             Order,
    const GeometryType&           rB0,
    const GeometryType&           rL0,
    const GeometryType&           rL1,
    const GeometryType&           rB1,
    const array_1d<double,3>&     rP00,
    const array_1d<double,3>&     rP01,
    const array_1d<double,3>&     rP10,
    const array_1d<double,3>&     rP11) const
{
    IntegrationPointsArrayType gp_list;
    gp_list.reserve(Order * Order);

    // 1-D Gauss nodes / weights on [0,1]
    std::vector<double> xi, w;
    GaussLegendreOnUnitInterval(Order, xi, w);

    for (std::size_t i=0; i<Order; ++i)
    for (std::size_t j=0; j<Order; ++j)
    {
        const double xi_i  = xi[i];
        const double eta_j = xi[j];          // same grid for η
        const double w_ij  = w[i] * w[j];

        // --- Jacobian ------------------------------------------------
        const array_1d<double,3> dXi =
            CoonsDerivativeFD(xi_i, eta_j, true ,
                              rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);

        const array_1d<double,3> dEta =
            CoonsDerivativeFD(xi_i, eta_j, false,
                              rB0,rL0,rL1,rB1, rP00,rP01,rP10,rP11);

        const array_1d<double,3> cross =
            MathUtils<double>::CrossProduct(dXi, dEta);

        const double jac = norm_2(cross);    // |X_ξ × X_η|

        // --- Global coordinates --------------------------------------
        const array_1d<double,3> X =
            CoonsPoint(xi_i, eta_j,
                       rB0,rL0,rL1,rB1,
                       rP00,rP01,rP10,rP11);

        // store *global* coords + weight (w·|J|)
        gp_list.emplace_back( IntegrationPoint<3>( X[0], X[1], X[2],
                                                   w_ij * jac ) );
    }

    return gp_list;
}


void SnakeExtendedSbmProcess::BuildParabolicNurbsData(
    Node::Pointer         pNode0,
    Node::Pointer         pNodeM,
    Node::Pointer         pNode2,
    PointerVector<Node>&  rCtrlPtsPointerVector,
    Vector&               rKnots,
    Vector&               rWeights,
    const double          tm)
{
    KRATOS_ERROR_IF(tm <= 0.0 || tm >= 1.0)
        << "tm must be in (0,1). Given: " << tm << std::endl;

    // ----------------------------------------------------------------
    // Control node 0 and 2 are the given ones
    // Control node 1 is created on-the-fly (ID = 0 → virtual, change if needed)
    // ----------------------------------------------------------------
    const array_1d<double,3> P0 = pNode0->Coordinates();
    const array_1d<double,3> Pm = pNodeM->Coordinates();
    const array_1d<double,3> P2 = pNode2->Coordinates();

    const double t  = tm;
    const double denom = 2.0 * t * (1.0 - t);

    array_1d<double,3> P1 =
          ( Pm
          - std::pow(1.0 - t, 2) * P0
          - std::pow(    t, 2)   * P2 ) / denom;

    // ----------------------------------------------------------------
    // Create a new node for P1 (ID = 0  → not stored in any ModelPart)
    // If you need it in a ModelPart, use pModelPart->CreateNewNode(...)
    // ----------------------------------------------------------------
    Node::Pointer pNode1 = Node::Pointer(new Node(1, P1));

    // ----------------------------------------------------------------
    // Fill PointerVector<Node>  (size 3)
    // ----------------------------------------------------------------
    rCtrlPtsPointerVector.clear();
    rCtrlPtsPointerVector.reserve(3);
    rCtrlPtsPointerVector.push_back(pNode0);
    rCtrlPtsPointerVector.push_back(pNode1);
    rCtrlPtsPointerVector.push_back(pNode2);

    // ----------------------------------------------------------------
    // Knot vector [0 0 0 1 1 1]  (Vector of size 6)
    // ----------------------------------------------------------------
    if (rKnots.size() != 6) rKnots.resize(6);
    rKnots[0] = rKnots[1] = rKnots[2] = 0.0;
    rKnots[3] = rKnots[4] = rKnots[5] = 1.0;

    // ----------------------------------------------------------------
    // Weights  (all 1 for polynomial parabola)
    // ----------------------------------------------------------------
    if (rWeights.size() != 3) rWeights.resize(3);
    rWeights[0] = rWeights[1] = rWeights[2] = 1.0;
}


void SnakeExtendedSbmProcess::FindClosestTruePointToSurrogateVertexByNurbs()
{
    // // get the data
    // mEchoLevel = mThisParameters["echo_level"].GetInt();
    // std::string iga_model_part_name = mThisParameters["model_part_name"].GetString();
    // std::string skin_model_part_name = mThisParameters["skin_model_part_name"].GetString();
    // std::string skin_model_part_inner_initial_name = mThisParameters["skin_model_part_inner_initial_name"].GetString();
    
    // // Loop over the surrogate vertex and call the search in radius or something like that
    // // TODO: Outer boyndary
    // const bool is_inner = true;
    
    // std::string surrogate_sub_model_part_name; 
    // std::string skin_sub_model_part_name; 
    // surrogate_sub_model_part_name = "surrogate_inner";
    // skin_sub_model_part_name = "inner";

    // ModelPart& r_extended_sbm_sub_model_part = mpIgaModelPart->CreateSubModelPart("extended_sbm");

    // ModelPart& r_skin_sub_model_part = mpSkinModelPart->GetSubModelPart(skin_sub_model_part_name);
    // ModelPart& r_surrogate_sub_model_part = mpIgaModelPart->GetSubModelPart(surrogate_sub_model_part_name);

    // auto p_surface = mpIgaModelPart->pGetGeometry(1);
    // IndexType id_brep_curve_on_surface = (mpIgaModelPart->GeometriesEnd()-1)->Id() + 1;
    
    // typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceType;
    // auto p_nurbs_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(
    //                         p_surface->pGetGeometryPart(Geometry<typename PointerVector<NodeType>::value_type>::BACKGROUND_GEOMETRY_INDEX));
    // IntegrationInfo surface_integration_info = p_nurbs_surface->GetDefaultIntegrationInfo();

    // // Get the mesh sizes from the surrogate model part
    // const Vector& knot_span_sizes = r_surrogate_sub_model_part.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);

    // double knot_span_reference_size = knot_span_sizes[0];
    // if (knot_span_sizes[1] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[1];}
    // if (knot_span_sizes.size() > 2) {if (knot_span_sizes[2] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[2];}}

    // std::string condition_name = "ExtendedSbmSolidCondition";
    // std::string element_name = "ExtendedSbmSolidElement";
    // // std::string element_name = "SolidElement";

    // // Loop over the nodes of the surrogate sub model part
    // IndexType iel = 1;
    // SizeType number_of_shape_functions_derivatives = 5;

    // IndexType first_condition_id = r_surrogate_sub_model_part.pGetElement(iel)->GetGeometry()[0].Id();
    // IndexType last_condition_id = r_surrogate_sub_model_part.pGetElement(iel)->GetGeometry()[1].Id();

    // SizeType size_surrogate_loop = last_condition_id - first_condition_id + 1;

    // for (SizeType j = 0; j < size_surrogate_loop; ++j) {
    //     auto p_brep_geometry = mpIgaModelPart->pGetGeometry(6 + j);
    //     auto p_brep_curve_on_surface_surrogate1_surrogate2 = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep_geometry);

    //     KRATOS_ERROR_IF(!p_brep_curve_on_surface_surrogate1_surrogate2) <<  ":::[SnakeExtendedSbmProcess]::: the geometry with id " << p_brep_curve_on_surface_surrogate1_surrogate2->Id() 
    //                                         << " is not a BrepCurveOnSurfaceType." << std::endl;

    //     NurbsInterval brep_domain_interval = p_brep_curve_on_surface_surrogate1_surrogate2->DomainInterval();
    //     CoordinatesArrayType surrogate_vertex_1 = ZeroVector(3); 
    //     CoordinatesArrayType surrogate_vertex_1_local_coords = ZeroVector(3);
    //     CoordinatesArrayType surrogate_vertex_2 = ZeroVector(3); 
    //     CoordinatesArrayType surrogate_vertex_2_local_coords = ZeroVector(3);
    //     surrogate_vertex_1_local_coords[0] = brep_domain_interval.GetT0();
    //     surrogate_vertex_2_local_coords[0] = brep_domain_interval.GetT1();

    //     p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_1, surrogate_vertex_1_local_coords);
    //     p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_2, surrogate_vertex_2_local_coords);

    //     // retrieve middle point of the brep
    //     CoordinatesArrayType surrogate_middle_point = ZeroVector(3); 
    //     CoordinatesArrayType surrogate_middle_point_local_coords = ZeroVector(3);
    //     surrogate_middle_point_local_coords[0] = 0.5 * (surrogate_vertex_1_local_coords[0] + surrogate_vertex_2_local_coords[0]);
    //     p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_middle_point, surrogate_middle_point_local_coords);

    //     IntegrationPoint<1> integration_point(surrogate_middle_point_local_coords[0]);
    //     IntegrationPointsArrayType surrogate_integration_points_list;
    //     surrogate_integration_points_list.push_back(integration_point);

    //     IntegrationInfo integration_info = p_brep_curve_on_surface_surrogate1_surrogate2->GetDefaultIntegrationInfo();
    //     GeometriesArrayType quadrature_point_list;
    //     p_brep_curve_on_surface_surrogate1_surrogate2->CreateQuadraturePointGeometries(quadrature_point_list, number_of_shape_functions_derivatives, 
    //                                                         surrogate_integration_points_list, integration_info);

    //     GeometryType::Pointer surrogate_brep_middle_geometry = quadrature_point_list(0);

    //     CoordinatesArrayType skin_vertex_1 = ZeroVector(3);
    //     CoordinatesArrayType skin_vertex_local_1 = ZeroVector(3);
    //     std::vector<array_1d<double, 3>> curve_derivatives;
    //     bool is_projected = this->ProjectToSkinBoundary(&*mpSkinModelPartInnerInitial, surrogate_vertex_1, skin_vertex_1, skin_vertex_local_1, curve_derivatives, 50);
        
    //     // second vertex
    //     CoordinatesArrayType skin_vertex_2 = ZeroVector(3);
    //     CoordinatesArrayType skin_vertex_local_2 = ZeroVector(3);
    //     std::vector<array_1d<double, 3>> curve_derivatives_2;
    //     bool is_projected_2 = this->ProjectToSkinBoundary(&*mpSkinModelPartInnerInitial, surrogate_vertex_2, skin_vertex_2, skin_vertex_local_2, curve_derivatives_2, 50);

    //     Vector active_range_knot_vector = ZeroVector(2);
    //     active_range_knot_vector[0] = 0;
    //     active_range_knot_vector[1] = 1;
    //     NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);

    //     // surrogate_1 - skin_1
    //     // create the brep connecting vertex and closest true point
    //     Point surrogate_1(surrogate_vertex_1);
    //     Point skin_1(skin_vertex_1);
        
    //     Node::Pointer p_surrogate1_brep_point = Node::Pointer(new Node(1, surrogate_1));
    //     Node::Pointer p_skin1_brep_point = Node::Pointer(new Node(2, skin_1));

    //     auto p_nurbs_curve_surrogate1_skin1 = this->CreateBrepCurve(p_surrogate1_brep_point, p_skin1_brep_point, active_range_knot_vector);
    //     auto p_brep_curve_surrogate1_skin1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_skin1);      
    
    //     // IntegrationInfo brep_integration_info_surrogate1_skin1 = p_brep_curve_surrogate1_skin1->GetDefaultIntegrationInfo();
    //     IntegrationPointsArrayType brep_integration_points_list_surrogate1_skin1;
    //     GeometriesArrayType brep_quadrature_point_list_surrogate1_skin1;

    //     p_brep_curve_surrogate1_skin1->CreateIntegrationPoints(brep_integration_points_list_surrogate1_skin1, integration_info);
    //     p_brep_curve_surrogate1_skin1->CreateQuadraturePointGeometries(brep_quadrature_point_list_surrogate1_skin1, number_of_shape_functions_derivatives, 
    //                                                         brep_integration_points_list_surrogate1_skin1, integration_info);

                                                            
        
    //     SizeType id = 1;
    //     if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
    //         id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;

    //     this->CreateConditions(
    //         brep_quadrature_point_list_surrogate1_skin1.ptr_begin(), brep_quadrature_point_list_surrogate1_skin1.ptr_end(),
    //         *mpIgaModelPart, condition_name, id, PropertiesPointerType(), knot_span_sizes, surrogate_brep_middle_geometry, "");
            

    //     // surrogate_2 - skin_2
    //     // create the brep connecting vertex and closest true point
    //     Point surrogate_2(surrogate_vertex_2);
    //     Point skin_2(skin_vertex_2);
        
    //     Node::Pointer p_surrogate2_brep_point = Node::Pointer(new Node(1, surrogate_2));
    //     Node::Pointer p_skin2_brep_point = Node::Pointer(new Node(2, skin_2));

    //     auto p_nurbs_curve_surrogate2_skin2 = this->CreateBrepCurve(p_surrogate2_brep_point, p_skin2_brep_point, active_range_knot_vector);
    //     auto p_brep_curve_surrogate2_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate2_skin2);      
    
    //     IntegrationPointsArrayType brep_integration_points_list_surrogate2_skin2;
    //     GeometriesArrayType brep_quadrature_point_list_surrogate2_skin2;

    //     p_brep_curve_surrogate2_skin2->CreateIntegrationPoints(brep_integration_points_list_surrogate2_skin2, integration_info);
    //     p_brep_curve_surrogate2_skin2->CreateQuadraturePointGeometries(brep_quadrature_point_list_surrogate2_skin2, number_of_shape_functions_derivatives, 
    //                                                 brep_integration_points_list_surrogate2_skin2, integration_info);
        
    //     if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
    //         id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;

    //     this->CreateConditions(
    //         brep_quadrature_point_list_surrogate2_skin2.ptr_begin(), brep_quadrature_point_list_surrogate2_skin2.ptr_end(),
    //         *mpIgaModelPart, condition_name, id, PropertiesPointerType(), knot_span_sizes, surrogate_brep_middle_geometry, "");


    //     // skin_1 - skin_2
    //     // create the brep connecting vertex and closest true point
    //     auto p_nurbs_curve_geometry = mpSkinModelPartInnerInitial->pGetGeometry(0);
    //     auto nurbs_curve_geometry = std::dynamic_pointer_cast<NurbsCurveGeometryType>(p_nurbs_curve_geometry);
    //     auto t0 = nurbs_curve_geometry->DomainInterval().GetT0();
    //     auto t1 = nurbs_curve_geometry->DomainInterval().GetT1();
    //     CoordinatesArrayType skin_middle_point_local_coords = ZeroVector(3);
    //     skin_middle_point_local_coords = (skin_vertex_local_1+skin_vertex_local_2)*0.5;

    //     if (skin_vertex_local_1[0] < skin_vertex_local_2[0]) 
    //         skin_middle_point_local_coords[0] = 0.5*((t1-t0) + skin_vertex_local_1[0]+skin_vertex_local_2[0]);

    //         if (skin_middle_point_local_coords[0] > t1) 
    //             skin_middle_point_local_coords[0] -= t1-t0;
    //         else if (skin_middle_point_local_coords[0] > r_skin_sub_model_part.Nodes().front().Id())
    //             skin_middle_point_local_coords[0] += t1-t0;

    //     CoordinatesArrayType skin_middle_point = ZeroVector(3);
    //     nurbs_curve_geometry->GlobalCoordinates(skin_middle_point, skin_middle_point_local_coords);
    //     Point p_skin_middle_point(skin_middle_point);
    //     Node::Pointer p_skin_middle_brep_point = Node::Pointer(new Node(2, p_skin_middle_point));

    //     PointerVector<Node> ctrl_pts;
    //     Vector               knots;
    //     Vector               weights;

    //     BuildParabolicNurbsData(p_skin1_brep_point, p_skin_middle_brep_point, p_skin2_brep_point,
    //                             ctrl_pts, knots, weights, 0.5);

    //     const double polynomial_degree = 2;
    //     typename NurbsCurveGeometry<3, PointerVector<Node>>::Pointer p_nurbs_curve_skin1_skin2(
    //             new NurbsCurveGeometry<3, PointerVector<Node>>(
    //                 ctrl_pts,
    //                 polynomial_degree,
    //                 knots,
    //                 weights)); 

    //     // auto p_nurbs_curve_skin1_skin2 = this->CreateBrepCurve(p_skin1_brep_point, p_skin2_brep_point, active_range_knot_vector);
    //     auto p_brep_curve_skin1_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_skin1_skin2);      
    
    //     IntegrationPointsArrayType brep_integration_points_list_skin1_skin2;
    //     GeometriesArrayType brep_quadrature_point_list_skin1_skin2;

    //     p_brep_curve_skin1_skin2->CreateIntegrationPoints(brep_integration_points_list_skin1_skin2, integration_info);

    //     const double p_brep_curve_skin1_skin2_length = (p_brep_curve_skin1_skin2->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX))->Length();
    //     for (auto& integration_point : brep_integration_points_list_skin1_skin2) {
    //         integration_point.SetWeight(integration_point.Weight() * p_brep_curve_skin1_skin2_length);
    //     }

    //     p_brep_curve_skin1_skin2->CreateQuadraturePointGeometries(brep_quadrature_point_list_skin1_skin2, number_of_shape_functions_derivatives, 
    //                                                               brep_integration_points_list_skin1_skin2, integration_info);
        
        
    //     if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
    //         id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;


    //     std::string identifier_name = "TRUE_BOUNDARY";
    //     this->CreateConditions(
    //         brep_quadrature_point_list_skin1_skin2.ptr_begin(), brep_quadrature_point_list_skin1_skin2.ptr_end(),
    //         r_extended_sbm_sub_model_part, condition_name, id, PropertiesPointerType(), knot_span_sizes, surrogate_brep_middle_geometry, identifier_name);

        
    //     // surrogate_1 - surrogate_2
    //     // create ONLY the brep connecting the two vertices
    //     auto p_nurbs_curve_surrogate1_surrogate2 = this->CreateBrepCurve(p_surrogate1_brep_point, p_surrogate2_brep_point, active_range_knot_vector);
    //     auto p_brep_curve_surrogate1_surrogate2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_surrogate2);      


        
    //     // Creation integration points on the cut elements

    //     IntegrationPointsArrayType surface_integration_points = CreateCoonsPatchGaussPoints(
    //                         3, /*Order*/
    //                         *p_brep_curve_surrogate1_surrogate2,   // B0
    //                         *p_brep_curve_surrogate1_skin1,       // L0
    //                         *p_brep_curve_surrogate2_skin2,       // L1
    //                         *p_brep_curve_skin1_skin2,            // B1
    //                         surrogate_1,  // P00
    //                         skin_1,       // P01
    //                         surrogate_2,  // P10
    //                         skin_2);      // P11
        
    //     GeometriesArrayType surface_quadrature_point_list;

    //     p_nurbs_surface->CreateQuadraturePointGeometries(surface_quadrature_point_list, number_of_shape_functions_derivatives, 
    //                                                     surface_integration_points, surface_integration_info);

    //     IndexType id_element = 1;
    //     if (mpIgaModelPart->GetRootModelPart().Elements().size() > 0)
    //         id_element = mpIgaModelPart->GetRootModelPart().Elements().back().Id() + 1;

    //     this->CreateElements(
    //         surface_quadrature_point_list.ptr_begin(), surface_quadrature_point_list.ptr_end(),
    //         *mpIgaModelPart, element_name, id_element, PropertiesPointerType(), surrogate_brep_middle_geometry);
    
    // }
    
}

}  // namespace Kratos.
