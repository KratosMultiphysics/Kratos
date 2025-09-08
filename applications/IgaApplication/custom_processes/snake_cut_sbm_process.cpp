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
#include "snake_cut_sbm_process.h"
#include "iga_application_variables.h"

namespace Kratos
{

SnakeCutSbmProcess::SnakeCutSbmProcess(
    Model& rModel, Parameters ThisParameters) : 
    SnakeSbmProcess(rModel, ThisParameters)
{

    KRATOS_ERROR_IF_NOT(ThisParameters.Has("cut_element_name")) << "::[SnakeCutSbmProcess]::" 
                    << "Missing \"cut_element_name\" section." << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("cut_interface_condition_name")) << "::[SnakeCutSbmProcess]::" 
                    << "Missing \"cut_interface_condition_name\" section." << std::endl;
    
    mpCutElementsSubModelPart = &(mpIgaModelPart->CreateSubModelPart("CutElements"));
    mpCutInterfaceSubModelPart = &(mpIgaModelPart->CreateSubModelPart("CutInterfaces"));
    mCutElementName = ThisParameters["cut_element_name"].GetString();
    mCutInterfaceConditionName = ThisParameters["cut_interface_condition_name"].GetString();
    mLambdaInner = 0.0;
    mLambdaOuter = 1.0;
    mInternalDivision = ThisParameters["number_internal_divisions"].GetInt();
}

void SnakeCutSbmProcess::CreateSbmExtendedGeometries()
{
    mEchoLevel = mThisParameters["echo_level"].GetInt();
    if (mpSkinModelPartInnerInitial->NumberOfNodes()>0 || mpSkinModelPartInnerInitial->NumberOfGeometries()>0) 
    {
        const auto& r_surrogate_sub_model_part_inner = mpIgaModelPart->GetSubModelPart("surrogate_inner");
        const auto& r_skin_sub_model_part_inner = mpSkinModelPart->GetSubModelPart("inner");

        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Creating the extended SBM geometries for the inner skin." << std::endl;
        CreateSbmExtendedGeometries<true>(r_skin_sub_model_part_inner, r_surrogate_sub_model_part_inner);
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Finished creating the extended SBM geometries for the inner skin." << std::endl;
    }
    if (mpSkinModelPartOuterInitial->NumberOfNodes()>0 || mpSkinModelPartOuterInitial->NumberOfGeometries()>0) 
    {
        const auto& r_surrogate_sub_model_part_outer = mpIgaModelPart->GetSubModelPart("surrogate_outer");
        const auto& r_skin_sub_model_part_outer = mpSkinModelPart->GetSubModelPart("outer");
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Creating the extended SBM geometries for the outer skin." << std::endl;
        CreateSbmExtendedGeometries<false>(r_skin_sub_model_part_outer, r_surrogate_sub_model_part_outer);
        KRATOS_INFO_IF("CreateSbmExtendedGeometries", mEchoLevel > 0)
            << "Finished creating the extended SBM geometries for the outer skin." << std::endl;
    }
}

template <bool TIsInnerLoop>
void SnakeCutSbmProcess::CreateSbmExtendedGeometries(
    const ModelPart& rSkinSubModelPart,
    const ModelPart& rSurrogateSubModelPart)
{
    // get the data

    // Create the testBins for the search in radius
    PointVector points;
    for (auto &i_node : rSkinSubModelPart.Nodes()) {
        points.push_back(PointTypePointer(new PointType(i_node.Id(), i_node.X(), i_node.Y(), i_node.Z())));
    }
    // Get the mesh sizes from the surrogate model part
    const Vector& knot_span_sizes = rSurrogateSubModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);

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

    //FIXME: create the test Bins for the interface nodes
    // Create the testBins for the search in radius
    PointVector points_interface;
    for (auto &i_node : rSkinSubModelPart.GetSubModelPart("interface_vertices").Nodes()) {
        points_interface.push_back(PointTypePointer(new PointType(i_node.Id(), i_node.X(), i_node.Y(), i_node.Z())));
    }

    DynamicBins testBins_interface(points_interface.begin(), points_interface.end());

    // Maximum number of results to be found in the search in radius
    const int number_of_results = 1e6; 
    const int number_of_results_interface = 10; 

    ModelPart::NodesContainerType::ContainerType results(number_of_results);
    std::vector<double> list_of_distances(number_of_results);

    BinSearchParameters bin_search_parameters(
        testBins, 
        number_of_results, 
        results, 
        list_of_distances, 
        search_radius);

    
    BinSearchParameters bin_search_interface_parameters(
        testBins_interface, 
        number_of_results_interface, 
        results, 
        list_of_distances, 
        search_radius);

    auto p_surface = mpIgaModelPart->pGetGeometry(1);
    IndexType id_brep_curve_on_surface = (mpIgaModelPart->GeometriesEnd()-1)->Id() + 1;
    
    auto p_nurbs_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(
                            p_surface->pGetGeometryPart(Geometry<typename PointerVector<NodeType>::value_type>::BACKGROUND_GEOMETRY_INDEX));
    IntegrationInfo surface_integration_info = p_nurbs_surface->GetDefaultIntegrationInfo();

    // FIXME: search the projection for each node of the surrogate model part
    SetSurrogateToSkinProjections<TIsInnerLoop>(rSurrogateSubModelPart, rSkinSubModelPart, bin_search_parameters, bin_search_interface_parameters);

    // Loop over the nodes of the surrogate sub model part
    IndexType iel = 1;
    SizeType number_of_shape_functions_derivatives = 2*p_nurbs_surface->PolynomialDegree(0)+1;

    IndexType first_condition_id;
    IndexType last_condition_id;
    IndexType starting_brep_id;
    SizeType size_surrogate_loop;

    const bool is_inner = TIsInnerLoop;
    if constexpr (TIsInnerLoop)  {
        first_condition_id = rSurrogateSubModelPart.pGetElement(iel)->GetGeometry()[0].Id();
        last_condition_id = rSurrogateSubModelPart.pGetElement(iel)->GetGeometry()[1].Id();
        size_surrogate_loop = last_condition_id - first_condition_id + 1;
        if (mpSkinModelPartOuterInitial->NumberOfNodes()>0 || mpSkinModelPartOuterInitial->NumberOfGeometries()>0) 
        {
            starting_brep_id = 2 + mpIgaModelPart->GetSubModelPart("surrogate_outer").NumberOfConditions(); //1 surface + outer surrogate loop
        }
        else
            starting_brep_id = 6; //1 surface + 4 external boundaries
    }
    else {
        size_surrogate_loop = rSurrogateSubModelPart.NumberOfConditions();
        first_condition_id = rSurrogateSubModelPart.ConditionsBegin()->Id();
        last_condition_id = first_condition_id + size_surrogate_loop - 1;
        starting_brep_id = 2; //1 surface 
    }

    for (SizeType j = 0; j < size_surrogate_loop; ++j) {
        auto p_brep_geometry = mpIgaModelPart->pGetGeometry(starting_brep_id + j);
        auto p_brep_curve_on_surface_surrogate1_surrogate2 = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep_geometry);

        KRATOS_ERROR_IF(!p_brep_curve_on_surface_surrogate1_surrogate2) <<  ":::[SnakeCutSbmProcess]::: the geometry with id " << p_brep_curve_on_surface_surrogate1_surrogate2->Id() 
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

        IntegrationParameters integration_parameters(
            number_of_shape_functions_derivatives, 
            integration_info, 
            knot_span_sizes);

        GeometriesArrayType quadrature_point_list;
        p_brep_curve_on_surface_surrogate1_surrogate2->CreateQuadraturePointGeometries(quadrature_point_list, number_of_shape_functions_derivatives, 
                                                            surrogate_integration_points_list, integration_info);

        GeometryType::Pointer surrogate_brep_middle_geometry = quadrature_point_list(0);

        // STORE THE SURROGATE BREP MIDDLE GEOMETRY FOR THE LATERAL BREPS
        const double tol = 1.0e-12;
        bool check_cond_1 = false;
        bool check_cond_2 = false;
        Node::Pointer p_surrogate_node_1;
        Node::Pointer p_surrogate_node_2;
        for (auto& r_surrogate_node : rSurrogateSubModelPart.Nodes())
        {
            // Check if this node coincides with surrogate_vertex_1
            if (norm_2(r_surrogate_node.Coordinates() - surrogate_vertex_1) < tol)
            {
                // Direct reference to the stored neighbour geometries
                auto& r_neighbour_geometries =
                    r_surrogate_node.GetValue(NEIGHBOUR_GEOMETRIES);
                
                p_surrogate_node_1 = &r_surrogate_node;

                // Append the new middle geometry; no SetValue call required
                r_neighbour_geometries.push_back(surrogate_brep_middle_geometry);

                check_cond_1 = true; // Only one node should satisfy the condition
            }
            else if (norm_2(r_surrogate_node.Coordinates() - surrogate_vertex_2) < tol)
            {
                // Direct reference to the stored neighbour geometries
                auto& r_neighbour_geometries =
                    r_surrogate_node.GetValue(NEIGHBOUR_GEOMETRIES);

                p_surrogate_node_2 = &r_surrogate_node;

                // Append the new middle geometry; no SetValue call required
                r_neighbour_geometries.push_back(surrogate_brep_middle_geometry);

                check_cond_2 = true; // Only one node should satisfy the condition
            }
            if (check_cond_1 && check_cond_2) {
                // If both conditions are satisfied, break the loop
                break;
            }
        }

        const double t0 = brep_domain_interval.GetT0();
        const double t1 = brep_domain_interval.GetT1();
        const double dt = (t1 - t0) / static_cast<double>(mInternalDivision);

        //------------------------------------------------------------------
        // 3. Loop over sub-intervals of the upper curve
        //------------------------------------------------------------------
        Node::Pointer p_first_node = p_surrogate_node_1;
        Node::Pointer p_second_node = nullptr;
        for (SizeType d = 0; d < mInternalDivision; ++d)
        {
            //---------------- 3.1 Current sub-interval --------------------
            const double sub_t0 = t0 +  d      * dt;
            const double sub_t1 = t0 + (d + 1) * dt;
            NurbsInterval sub_interval(sub_t0, sub_t1);

            surrogate_vertex_1_local_coords[0] = sub_t0;
            surrogate_vertex_2_local_coords[0] = sub_t1;

            p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_1, surrogate_vertex_1_local_coords);
            p_brep_curve_on_surface_surrogate1_surrogate2->GlobalCoordinates(surrogate_vertex_2, surrogate_vertex_2_local_coords);

            if (d != mInternalDivision - 1) {
                // create the new node 
                p_second_node = Node::Pointer(new Node(0, surrogate_vertex_2));
            }
            else {
                // last division -> use the second surrogate node
                p_second_node = p_surrogate_node_2;
            }

            CreateCutAndSkinQuadraturePoints<TIsInnerLoop>(integration_parameters, bin_search_parameters,
                                        p_nurbs_surface, p_first_node, p_second_node,
                                        surrogate_brep_middle_geometry,
                                        *mpIgaModelPart, rSkinSubModelPart);    

            p_first_node = p_second_node;
        }
    }
    
    //---------------------------------------------------------------------------
    bool is_entering = false;
    for (IndexType i_cond_id = first_condition_id; i_cond_id <= last_condition_id; ++i_cond_id)
    {
        is_entering = !is_entering;

        const auto& surrogate_condition = rSurrogateSubModelPart.pGetCondition(i_cond_id);

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

            const IndexType id_closest_true_node = p_surrogate_1->GetValue(PROJECTION_NODE_ID);
            
            const auto skin_vertex_1 = rSkinSubModelPart.pGetNode(id_closest_true_node);

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
                *mpCutInterfaceSubModelPart, mCutInterfaceConditionName, id, PropertiesPointerType(), knot_span_sizes, neighbour_geometries);

        }

        if (!is_second_surrogate_already_computed) {
            p_surrogate_2->SetValue(ACTIVATION_LEVEL, 1.0); 
           
            const IndexType id_closest_true_node = p_surrogate_2->GetValue(PROJECTION_NODE_ID);
            const auto skin_vertex_2 = rSkinSubModelPart.pGetNode(id_closest_true_node);


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
                *mpCutInterfaceSubModelPart, mCutInterfaceConditionName, id, PropertiesPointerType(), knot_span_sizes, neighbour_geometries);

        }
    }
}

template <bool TIsInnerLoop>
void SnakeCutSbmProcess::CreateCutAndSkinQuadraturePoints(
    IntegrationParameters& rIntegrationParameters,
    BinSearchParameters& rBinSearchParameters,
    NurbsSurfaceType::Pointer& pNurbsSurface,
    const Node::Pointer& pSurrogateNode1, 
    const Node::Pointer& pSurrogateNode2, 
    const GeometryType::Pointer& rSurrogateBrepMiddleGeometry,
    ModelPart& rIgaModelPart,
    const ModelPart& rSkinSubModelPart)
{
    const IndexType id_closest_true_node = pSurrogateNode1->GetValue(PROJECTION_NODE_ID);

    const auto p_skin_node_1 = rSkinSubModelPart.pGetNode(id_closest_true_node);
    
    // search the projection of the second vertex
    DynamicBinsPointerType p_surrogate_vertex_2 = DynamicBinsPointerType(new PointType(1, pSurrogateNode2->X(), pSurrogateNode2->Y(), pSurrogateNode2->Z()));

    const IndexType id_closest_true_node_2 = pSurrogateNode2->GetValue(PROJECTION_NODE_ID);

    const auto p_skin_node_2 = rSkinSubModelPart.pGetNode(id_closest_true_node_2);

    // retrieve condition name for the skin condition
    auto connected_layers_1 = pSurrogateNode1->GetValue(CONNECTED_LAYERS);
    auto connected_layers_2 = pSurrogateNode2->GetValue(CONNECTED_LAYERS);
    std::string layer_name = "";   
    std::string condition_name = "";
    IndexType condition_count = 0;
    bool layer_found = false;
    // Find the common layer between the two surrogate nodes
    for (auto& layer : connected_layers_1) {
        for (auto& layer_2 : connected_layers_2) {
            if (layer == layer_2) {
                layer_name = layer;
                condition_name = pSurrogateNode1->GetValue(CONNECTED_CONDITIONS)[condition_count];
                layer_found = true;
                break;                  
            }
        }
        condition_count++;
    }

    KRATOS_ERROR_IF(!layer_found) << ":::[SnakeCutSbmProcess]::: No common layer found between the two surrogate nodes "
                                    << pSurrogateNode1->Id() << " and " << pSurrogateNode2->Id() <<  "\n"
                                    << pSurrogateNode1->Coordinates() << pSurrogateNode2->Coordinates() << "\n"
                                    << connected_layers_1 << connected_layers_2 << std::endl;

    ModelPart& r_layer_model_part = rIgaModelPart.HasSubModelPart(layer_name) ? 
                                    rIgaModelPart.GetSubModelPart(layer_name) : 
                                    rIgaModelPart.CreateSubModelPart(layer_name);       
                                    

    Vector active_range_knot_vector = ZeroVector(2);
    active_range_knot_vector[0] = 0;
    active_range_knot_vector[1] = 1;
    NurbsInterval brep_active_range(active_range_knot_vector[0], active_range_knot_vector[1]);

    // surrogate_1 - skin_1
    // create the brep connecting vertex and closest true point
    Point skin_1(*p_skin_node_1);
    
    Node::Pointer p_skin1_brep_point = Node::Pointer(new Node(2, skin_1));
    
    auto p_nurbs_curve_surrogate1_skin1 = this->CreateBrepCurve(pSurrogateNode1, p_skin1_brep_point, active_range_knot_vector);
    auto p_brep_curve_surrogate1_skin1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_skin1);      
                                            
    // surrogate_2 - skin_2
    // create the brep connecting vertex and closest true point
    Point skin_2(*p_skin_node_2);
    
    Node::Pointer p_skin2_brep_point = Node::Pointer(new Node(2, skin_2));

    auto p_nurbs_curve_surrogate2_skin2 = this->CreateBrepCurve(pSurrogateNode2, p_skin2_brep_point, active_range_knot_vector);
    auto p_brep_curve_surrogate2_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate2_skin2);      

    // skin_1 - skin_2
    //FIXME:
    // auto p_nurbs_curve_skin1_skin2 = FitQuadraticBezierBetween<TIsInnerLoop>(
    //     rSkinSubModelPart, id_closest_true_condition, id_closest_true_condition_2,
    //     p_skin1_brep_point, p_skin2_brep_point);
    

    const int p = 2;
    auto p_nurbs_curve_skin1_skin2 = this->CreateBrepCurve(p_skin1_brep_point, p_skin2_brep_point, active_range_knot_vector);

    // if (norm_2(skin_2 - skin_1) > rBinSearchParameters.SearchRadius/10.0)
        // p_nurbs_curve_skin1_skin2 = FitUV_BetweenSkinConditions_Generic<TIsInnerLoop>(
        //     rSkinSubModelPart, *pNurbsSurface, id_closest_true_condition, id_closest_true_condition_2, p, /*ridge=*/1e-14);


    auto p_brep_curve_skin1_skin2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_skin1_skin2);      

    IntegrationPointsArrayType brep_integration_points_list_skin1_skin2;
    GeometriesArrayType brep_quadrature_point_list_skin1_skin2;

    p_brep_curve_skin1_skin2->CreateIntegrationPoints(brep_integration_points_list_skin1_skin2, rIntegrationParameters.CurveIntegrationInfo);

    const double p_brep_curve_skin1_skin2_length = (p_brep_curve_skin1_skin2->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX))->Length();
    for (auto& integration_point : brep_integration_points_list_skin1_skin2) {
        integration_point.SetWeight(integration_point.Weight() * p_brep_curve_skin1_skin2_length);
    }

    p_brep_curve_skin1_skin2->CreateQuadraturePointGeometries(brep_quadrature_point_list_skin1_skin2, rIntegrationParameters.NumberOfShapeFunctionsDerivatives, 
                                                            brep_integration_points_list_skin1_skin2, rIntegrationParameters.CurveIntegrationInfo);
    
    SizeType id = 1;
    if (mpIgaModelPart->GetRootModelPart().Conditions().size() > 0)
        id = mpIgaModelPart->GetRootModelPart().Conditions().back().Id() + 1;
    
    std::vector<Geometry<Node>::Pointer> neighbour_geometries_skin1_skin2;
    neighbour_geometries_skin1_skin2.push_back(rSurrogateBrepMiddleGeometry);
    
    this->CreateConditions(
        brep_quadrature_point_list_skin1_skin2.ptr_begin(), brep_quadrature_point_list_skin1_skin2.ptr_end(),
        r_layer_model_part, condition_name, id, PropertiesPointerType(), rIntegrationParameters.KnotSpanSizes, neighbour_geometries_skin1_skin2);

    // FIXME: check for void elements/true coincident with surrogate boundary
    if (
        (std::abs(pSurrogateNode1->X() - p_skin_node_1->X()) < 1e-14 &&
         std::abs(pSurrogateNode1->X() - p_skin_node_2->X()) < 1e-14 &&
         std::abs(pSurrogateNode1->X() - pSurrogateNode2->X()) < 1e-14)  ||
        (std::abs(pSurrogateNode1->Y() - p_skin_node_1->Y()) < 1e-14 &&
         std::abs(pSurrogateNode1->Y() - p_skin_node_2->Y()) < 1e-14 &&
         std::abs(pSurrogateNode1->Y() - pSurrogateNode2->Y()) < 1e-14))
        return;
    
    // surrogate_1 - surrogate_2
    // create ONLY the brep connecting the two vertices
    auto p_nurbs_curve_surrogate1_surrogate2 = this->CreateBrepCurve(pSurrogateNode1, pSurrogateNode2, active_range_knot_vector);
    auto p_brep_curve_surrogate1_surrogate2 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate1_surrogate2);      


    // Creation integration points on the cut elements
    GeometriesArrayType surface_quadrature_point_list;
    IntegrationInfo surface_integration_info = pNurbsSurface->GetDefaultIntegrationInfo();
    if constexpr (TIsInnerLoop)  
    {

        // Bottom B0: surrogate_1 — surrogate_2  (you already have the Brep curve)
        auto p_uv_B0 = MakeUV_From3DEdge(*p_brep_curve_surrogate1_surrogate2, *pNurbsSurface /*, n_samples=9..13*/);

        // Left  L0: surrogate_1 — skin_1
        auto p_uv_L0 = MakeUV_From3DEdge(*p_brep_curve_surrogate1_skin1, *pNurbsSurface);

        // Right L1: surrogate_2 — skin_2
        auto p_uv_L1 = MakeUV_From3DEdge(*p_brep_curve_surrogate2_skin2, *pNurbsSurface);

        // auto p_uv_B1 = FitUV_BetweenSkinConditions_Generic<TIsInnerLoop>(
        //     rSkinSubModelPart, *pNurbsSurface, id_closest_true_condition, id_closest_true_condition_2, p, /*ridge=*/1e-14);

        // auto p_uv_B1 = FitQuadraticBezierBetween<TIsInnerLoop>(rSkinSubModelPart,
        //     id_closest_true_condition, id_closest_true_condition_2, p_skin1_brep_point, p_skin2_brep_poin

        auto p_uv_B1 = MakeUV_From3DEdge(*p_brep_curve_skin1_skin2, *pNurbsSurface);

        // 2) UV-Coons Gauss points (increase Order as needed, e.g. 5–7)
        IntegrationPointsArrayType surface_integration_points =
            CreateCoonsPatchGaussPointsUV(/*Order=*/4, *p_uv_B0, *p_uv_L0, *p_uv_L1, *p_uv_B1);

        pNurbsSurface->CreateQuadraturePointGeometries(surface_quadrature_point_list, rIntegrationParameters.NumberOfShapeFunctionsDerivatives, 
                                                        surface_integration_points, surface_integration_info);
    }
    else
    {
        auto p_nurbs_curve_skin2_skin1 = ReverseBezierUV_Generic(p_nurbs_curve_skin1_skin2);
        auto p_brep_curve_skin2_skin1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_skin2_skin1);  
        
        // // // surrogate_2 - surrogate_1
        auto p_nurbs_curve_surrogate2_surrogate1 = this->CreateBrepCurve(pSurrogateNode2, pSurrogateNode1, active_range_knot_vector);
        auto p_brep_curve_surrogate2_surrogate1 = Kratos::make_shared<BrepCurveType>(p_nurbs_curve_surrogate2_surrogate1);   

        auto p_uv_B0 = MakeUV_From3DEdge(*p_brep_curve_surrogate2_surrogate1, *pNurbsSurface /*, n_samples=9..13*/);

        // Left  L0: surrogate_1 — skin_1
        // auto p_brep_curve_surrogate1_skin1 = /* as in your code */;
        auto p_uv_L0 = MakeUV_From3DEdge(*p_brep_curve_surrogate1_skin1, *pNurbsSurface);

        // Right L1: surrogate_2 — skin_2
        // auto p_brep_curve_surrogate2_skin2 = /* as in your code */;
        auto p_uv_L1 = MakeUV_From3DEdge(*p_brep_curve_surrogate2_skin2, *pNurbsSurface);

        // // Top   B1: skin_1 — skin_2  (use all real skin points between the two conditions)
        // auto p_uv_B1 = FitUV_BetweenSkinConditions<TIsInnerLoop>(rSkinSubModelPart, *pNurbsSurface,
        //     id_closest_true_condition, id_closest_true_condition_2);

        auto p_uv_B1 = MakeUV_From3DEdge(*p_brep_curve_skin2_skin1, *pNurbsSurface);

        // 2) UV-Coons Gauss points (increase Order as needed, e.g. 5–7)
        IntegrationPointsArrayType surface_integration_points =
            CreateCoonsPatchGaussPointsUV(/*Order=*/4, *p_uv_B0, *p_uv_L1, *p_uv_L0, *p_uv_B1);
                                                        

        pNurbsSurface->CreateQuadraturePointGeometries(surface_quadrature_point_list, rIntegrationParameters.NumberOfShapeFunctionsDerivatives, 
                                                        surface_integration_points, surface_integration_info);
    }

    IndexType id_element = 1;
    if (mpCutElementsSubModelPart->GetRootModelPart().Elements().size() > 0)
        id_element = mpCutElementsSubModelPart->GetRootModelPart().Elements().back().Id() + 1;

    this->CreateElements(
        surface_quadrature_point_list.ptr_begin(), surface_quadrature_point_list.ptr_end(),
        *mpCutElementsSubModelPart, mCutElementName, id_element, PropertiesPointerType(), neighbour_geometries_skin1_skin2);
}






bool SnakeCutSbmProcess::ProjectToSkinBoundary(
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


void SnakeCutSbmProcess::CreateConditions(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    const std::string& rConditionName,
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


void SnakeCutSbmProcess::CreateElements(
    typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
    typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
    ModelPart& rModelPart,
    const std::string& rElementName,
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
void SnakeCutSbmProcess::GaussLegendreOnUnitInterval(
    const std::size_t      Order,
    std::vector<double>&   rXi,
    std::vector<double>&   rWeight)
{
    KRATOS_ERROR_IF(Order < 1 || Order > 10)
        << "Gauss order " << Order << " not implemented (1…10 supported)." << std::endl;

    // Gauss–Legendre nodes/weights on [-1,1] (rows: n=1..10).
    // Only the first 'n' entries of each row are used; the rest are zeros.
    static const double sXi[10][10] = {
        // n = 1
        {  0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        // n = 2
        { -0.5773502691896257,  +0.5773502691896257, 0, 0, 0, 0, 0, 0, 0, 0 },
        // n = 3
        {  0.0, -0.7745966692414834, +0.7745966692414834, 0, 0, 0, 0, 0, 0, 0 },
        // n = 4
        { -0.3399810435848563, +0.3399810435848563, -0.8611363115940526, +0.8611363115940526, 0, 0, 0, 0, 0, 0 },
        // n = 5
        {  0.0, -0.5384693101056831, +0.5384693101056831, -0.9061798459386640, +0.9061798459386640, 0, 0, 0, 0, 0 },
        // n = 6
        { -0.2386191860831969, +0.2386191860831969, -0.6612093864662645, +0.6612093864662645,
          -0.9324695142031521, +0.9324695142031521, 0, 0, 0, 0 },
        // n = 7
        {  0.0, -0.4058451513773972, +0.4058451513773972, -0.7415311855993945, +0.7415311855993945,
          -0.9491079123427585, +0.9491079123427585, 0, 0, 0 },
        // n = 8
        { -0.1834346424956498, +0.1834346424956498, -0.5255324099163290, +0.5255324099163290,
          -0.7966664774136267, +0.7966664774136267, -0.9602898564975363, +0.9602898564975363, 0, 0 },
        // n = 9
        {  0.0, -0.3242534234038089, +0.3242534234038089, -0.6133714327005904, +0.6133714327005904,
          -0.8360311073266358, +0.8360311073266358, -0.9681602395076261, +0.9681602395076261, 0 },
        // n = 10
        { -0.1488743389816312, +0.1488743389816312, -0.4333953941292472, +0.4333953941292472,
          -0.6794095682990244, +0.6794095682990244, -0.8650633666889845, +0.8650633666889845,
          -0.9739065285171717, +0.9739065285171717 }
    };

    static const double sW[10][10] = {
        // n = 1
        {  2.0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        // n = 2
        {  1.0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0 },
        // n = 3
        {  0.8888888888888888, 0.5555555555555556, 0.5555555555555556, 0, 0, 0, 0, 0, 0, 0 },
        // n = 4
        {  0.6521451548625461, 0.6521451548625461, 0.3478548451374539, 0.3478548451374539, 0, 0, 0, 0, 0, 0 },
        // n = 5
        {  0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891, 0, 0, 0, 0, 0 },
        // n = 6
        {  0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.3607615730481386,
           0.1713244923791704, 0.1713244923791704, 0, 0, 0, 0 },
        // n = 7
        {  0.4179591836734694, 0.3818300505051189, 0.3818300505051189, 0.2797053914892766, 0.2797053914892766,
           0.1294849661688697, 0.1294849661688697, 0, 0, 0 },
        // n = 8
        {  0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.3137066458778873,
           0.2223810344533745, 0.2223810344533745, 0.1012285362903763, 0.1012285362903763, 0, 0 },
        // n = 9
        {  0.3302393550012598, 0.3123470770400029, 0.3123470770400029, 0.2606106964029354, 0.2606106964029354,
           0.1806481606948574, 0.1806481606948574, 0.0812743883615744, 0.0812743883615744, 0 },
        // n = 10
        {  0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2692667193099963,
           0.2190863625159820, 0.2190863625159820, 0.1494513491505806, 0.1494513491505806,
           0.0666713443086881, 0.0666713443086881 }
    };

    const std::size_t n = Order;

    rXi.resize(n);
    rWeight.resize(n);

    for (std::size_t k = 0; k < n; ++k) {
        // Map node from [-1,1] to [0,1]
        rXi[k]     = 0.5 * (sXi[n-1][k] + 1.0);
        // Rescale weight: w' = w / 2
        rWeight[k] = 0.5 *  sW [n-1][k];
    }
}
// ------------------------------------------------------------------
// Global point on Brep curve
// ------------------------------------------------------------------
array_1d<double,3> SnakeCutSbmProcess::GlobalPoint(
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
array_1d<double,3> SnakeCutSbmProcess::CoonsPoint(
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
array_1d<double,3> SnakeCutSbmProcess::CoonsDerivativeFD(
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
SnakeCutSbmProcess::IntegrationPointsArrayType
SnakeCutSbmProcess::CreateCoonsPatchGaussPoints(
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


void SnakeCutSbmProcess::BuildParabolicNurbsData(
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

template <bool TIsInnerLoop>
void SnakeCutSbmProcess::SetSurrogateToSkinProjections(
    const ModelPart& rSurrogateSubModelPart, 
    const ModelPart& rSkinSubModelPart,
    BinSearchParameters& rBinSearchParameters, 
    BinSearchParameters& rBinSearchInterfaceParameters)
{
    auto has_projection = [](Node& rNode){
        return rNode.Has(PROJECTION_NODE_ID);
    };

    bool is_entering = false;

    for (auto& r_surrogate_condition : rSurrogateSubModelPart.Conditions())
    {
        // Search the closest condition in the skin model part
        is_entering = !is_entering;

        Node::Pointer p_surrogate_node_1 = r_surrogate_condition.pGetGeometry()->pGetPoint(0);
        Node::Pointer p_surrogate_node_2 = r_surrogate_condition.pGetGeometry()->pGetPoint(1);

        const IndexType first_node_id = rSkinSubModelPart.NodesBegin()->Id();
        const IndexType last_node_id  = first_node_id + rSkinSubModelPart.NumberOfNodes() - 1;
        const IndexType size_node_loop = rSkinSubModelPart.NumberOfNodes();
        
        const bool has_proj_1 = has_projection(*p_surrogate_node_1);
        const bool has_proj_2 = has_projection(*p_surrogate_node_2);

        IndexType skin_node_id_1 = -1;
        IndexType skin_node_id_2 = -1;

        if (has_proj_1 && has_proj_2) {
            auto node_1_connected_layers = p_surrogate_node_1->GetValue(CONNECTED_LAYERS);
            auto node_2_connected_layers = p_surrogate_node_2->GetValue(CONNECTED_LAYERS);

            bool have_common_layer = false;
            for (const auto& layer1 : node_1_connected_layers) {
                if (std::find(node_2_connected_layers.begin(),
                            node_2_connected_layers.end(),
                            layer1) != node_2_connected_layers.end()) {
                    have_common_layer = true;
                    break;
                }
            }
            if (have_common_layer) {
                // go to the next surrogate condition
                if (p_surrogate_node_1->GetValue(PROJECTION_NODE_ID) == p_surrogate_node_2->GetValue(PROJECTION_NODE_ID)) {
                    IndexType id_node = p_surrogate_node_1->GetValue(PROJECTION_NODE_ID);
                    if constexpr (TIsInnerLoop) {
                        const IndexType wrap_prev =
                            ( (id_node - 1 - first_node_id + size_node_loop) % size_node_loop ) + first_node_id;
                        p_surrogate_node_2->SetValue(PROJECTION_NODE_ID, wrap_prev);
                    } else {
                        const IndexType wrap_next =
                            ( (id_node + 1 - first_node_id + size_node_loop) % size_node_loop ) + first_node_id;
                        p_surrogate_node_2->SetValue(PROJECTION_NODE_ID, wrap_next);
                    }
                } 
                continue;
            }
            else
            {
                AssestProjectionsFeasibility(rSkinSubModelPart, p_surrogate_node_1, p_surrogate_node_2, rBinSearchInterfaceParameters);
                continue;
            }
        }

        std::vector<std::string> forced_layers;

        if (has_proj_1 && !has_proj_2) {
            forced_layers = p_surrogate_node_1->GetValue(CONNECTED_LAYERS); //FIXME:
            skin_node_id_1 = p_surrogate_node_1->GetValue(PROJECTION_NODE_ID);
        } else if (!has_proj_1 && has_proj_2) {
            forced_layers = p_surrogate_node_2->GetValue(CONNECTED_LAYERS);
            skin_node_id_2 = p_surrogate_node_2->GetValue(PROJECTION_NODE_ID);
        }

        if (!has_proj_1)
        {
            DynamicBinsPointerType p_surrogate_vertex_1 = DynamicBinsPointerType(new PointType(1, p_surrogate_node_1->X(), p_surrogate_node_1->Y(), p_surrogate_node_1->Z()));
            
            rBinSearchParameters.reset();
            SizeType obtained_results = rBinSearchParameters.TestBins.SearchInRadius(
                                                        *p_surrogate_vertex_1, 
                                                        rBinSearchParameters.SearchRadius, 
                                                        rBinSearchParameters.Results.begin(), 
                                                        rBinSearchParameters.ListOfDistances.begin(), 
                                                        rBinSearchParameters.NumberOfResults);
            double minimum_distance = 1e14;
            // Find the nearest node
            IndexType nearest_node_id;
            for (IndexType k = 0; k < obtained_results; k++) {
                double current_distance = rBinSearchParameters.ListOfDistances[k];   
                if (current_distance < minimum_distance) { 
                    IndexType current_id = rBinSearchParameters.Results[k]->Id();
                    // check orientation
                    if (skin_node_id_2 != -1)
                    {                  
                        if (is_entering)
                        {
                            if constexpr (TIsInnerLoop) {
                                if ((current_id > skin_node_id_2 && double(abs(skin_node_id_2-current_id)) < size_node_loop/2)
                                    || current_id == skin_node_id_2) {
                                    continue;
                                }
                            }
                            else {
                                if ((current_id < skin_node_id_2 && double(abs(skin_node_id_2-current_id)) < size_node_loop/2 )
                                    || current_id == skin_node_id_2) {
                                    continue;
                                }
                            }
                        }
                        else 
                        {
                            if constexpr (TIsInnerLoop) {
                                if ((current_id < skin_node_id_2 && double(abs(skin_node_id_2-current_id)) < size_node_loop/2)
                                    || current_id == skin_node_id_2) {
                                    continue;
                                }
                            }
                            else {
                                if ((current_id > skin_node_id_2 && double(abs(skin_node_id_2-current_id)) < size_node_loop/2 )
                                    || current_id == skin_node_id_2) {
                                    continue;
                                }
                            }
                        }
                    }

                    minimum_distance = current_distance;
                    nearest_node_id = current_id;
                }
            }
            KRATOS_ERROR_IF(obtained_results == 0) << "::[SnakeSbmProcess]:: Zero points found in serch for projection of point: " <<
                p_surrogate_vertex_1 << std::endl;

            // set the values to the node
            skin_node_id_1 = nearest_node_id;

            p_surrogate_node_1->SetValue(PROJECTION_NODE_ID, skin_node_id_1);

            auto connected_layers = rSkinSubModelPart.GetNode(skin_node_id_1).GetValue(CONNECTED_LAYERS);
            auto connected_conditions = rSkinSubModelPart.GetNode(skin_node_id_1).GetValue(CONNECTED_CONDITIONS);

            p_surrogate_node_1->SetValue(CONNECTED_LAYERS, connected_layers);
            p_surrogate_node_1->SetValue(CONNECTED_CONDITIONS, connected_conditions);

            // check feasibility
            bool is_feasible = false;
            auto nearest_node_connected_layers = rSkinSubModelPart.GetNode(nearest_node_id).GetValue(CONNECTED_LAYERS);
            if (forced_layers.size()>0) {
                for (const auto& forced_layer : forced_layers) {
                    if (std::find(nearest_node_connected_layers.begin(), nearest_node_connected_layers.end(), forced_layer) != nearest_node_connected_layers.end()) {
                        is_feasible = true;
                        break;
                    }
                }
            } else {
                is_feasible = true;
            }
            
            if (!is_feasible) 
                AssestProjectionsFeasibility(rSkinSubModelPart, p_surrogate_node_1, p_surrogate_node_2, rBinSearchInterfaceParameters);

            skin_node_id_1 = p_surrogate_node_1->GetValue(PROJECTION_NODE_ID);
            forced_layers = rSkinSubModelPart.GetNode(skin_node_id_1).GetValue(CONNECTED_LAYERS);
        }
        
        if (!has_proj_2)
        {
            // search the projection of the second vertex
            DynamicBinsPointerType p_surrogate_vertex_2 = DynamicBinsPointerType(new PointType(1, p_surrogate_node_2->X(), p_surrogate_node_2->Y(), p_surrogate_node_2->Z()));
            
            // Apply the SearchInRadius
            rBinSearchParameters.reset();
            SizeType obtained_results = rBinSearchParameters.TestBins.SearchInRadius(
                                                        *p_surrogate_vertex_2, 
                                                        rBinSearchParameters.SearchRadius, 
                                                        rBinSearchParameters.Results.begin(), 
                                                        rBinSearchParameters.ListOfDistances.begin(), 
                                                        rBinSearchParameters.NumberOfResults);                                     
            double minimum_distance = 1e14;
            IndexType nearest_node_id = -1;

            KRATOS_ERROR_IF(obtained_results == 0) << ":::[SnakeCutSbmProcess]::: No condition found in the skin model part for the surrogate node "
                                            << p_surrogate_node_2->Id() << " at coordinates " << p_surrogate_node_2->Coordinates() << std::endl;

            for (IndexType k = 0; k < obtained_results; k++) {
                double current_distance = rBinSearchParameters.ListOfDistances[k];   
               
                if (current_distance < minimum_distance) { 

                    IndexType current_id = rBinSearchParameters.Results[k]->Id();
                    if (is_entering)
                    {
                        if constexpr (TIsInnerLoop) {
                            if ((current_id < skin_node_id_1 && double(abs(skin_node_id_1-current_id)) < size_node_loop/2)
                                || current_id == skin_node_id_1) {
                                continue;
                            }
                        }
                        else {
                            if ((current_id > skin_node_id_1 && double(abs(skin_node_id_1-current_id)) < size_node_loop/2 )
                                || current_id == skin_node_id_1) {
                                continue;
                            }
                        }
                    }
                    else 
                    {
                        if constexpr (TIsInnerLoop) {
                            if ((current_id > skin_node_id_1 && double(abs(skin_node_id_1-current_id)) < size_node_loop/2)
                                || current_id == skin_node_id_1) {
                                continue;
                            }
                        }
                        else {
                            if ((current_id < skin_node_id_1 && double(abs(skin_node_id_1-current_id)) < size_node_loop/2 )
                                || current_id == skin_node_id_1) {
                                continue;
                            }
                        }
                    }
                    minimum_distance = current_distance;
                    nearest_node_id = current_id;
                }
            }

            KRATOS_ERROR_IF(nearest_node_id == -1) << "::[SnakeSbmProcess]:: No suitable second node found in search for projection of point: " <<
                p_surrogate_vertex_2 << " with first node id: " << skin_node_id_1 << std::endl;

            KRATOS_ERROR_IF(obtained_results == 0) << "::[SnakeSbmProcess]:: Zero points found in serch for projection of point: " <<
                p_surrogate_vertex_2 << std::endl;

            // set the values to the node
            skin_node_id_2 = nearest_node_id;

            p_surrogate_node_2->SetValue(PROJECTION_NODE_ID, skin_node_id_2);

            auto connected_layers = rSkinSubModelPart.GetNode(skin_node_id_2).GetValue(CONNECTED_LAYERS);
            auto connected_conditions = rSkinSubModelPart.GetNode(skin_node_id_2).GetValue(CONNECTED_CONDITIONS);

            p_surrogate_node_2->SetValue(CONNECTED_LAYERS, connected_layers);
            p_surrogate_node_2->SetValue(CONNECTED_CONDITIONS, connected_conditions);

            // check feasibility
            bool is_feasible = false;
            auto nearest_node_connected_layers = rSkinSubModelPart.GetNode(nearest_node_id).GetValue(CONNECTED_LAYERS);
            if (forced_layers.size()>0) {
                for (const auto& forced_layer : forced_layers) {
                    if (std::find(nearest_node_connected_layers.begin(), nearest_node_connected_layers.end(), forced_layer) != nearest_node_connected_layers.end()) {
                        is_feasible = true;
                        break;
                    }
                }
            } else {
                is_feasible = true;
            }

            if (!is_feasible) 
                AssestProjectionsFeasibility(rSkinSubModelPart, p_surrogate_node_1, p_surrogate_node_2, rBinSearchInterfaceParameters);
        }
    }
}

void SnakeCutSbmProcess::AssestProjectionsFeasibility(
    const ModelPart& rSkinSubModelPart,
    Node::Pointer pSurrogateNode1, 
    Node::Pointer pSurrogateNode2,
    BinSearchParameters& rBinSearchInterfaceParameters)
{
    // fix the projections 
    // keep the first projection and change the second one
    // search the projection of the second vertex
    IndexType skin_node_id_1 = pSurrogateNode1->GetValue(PROJECTION_NODE_ID);
    std::vector<std::string> forced_layers = pSurrogateNode1->GetValue(CONNECTED_LAYERS);

    DynamicBinsPointerType p_surrogate_vertex_2 = DynamicBinsPointerType(new PointType(1, pSurrogateNode2->X(), pSurrogateNode2->Y(), pSurrogateNode2->Z()));
    
    // Apply the SearchInRadius
    rBinSearchInterfaceParameters.reset();
    SizeType obtained_results_2 = rBinSearchInterfaceParameters.TestBins.SearchInRadius(
                                                *p_surrogate_vertex_2, 
                                                rBinSearchInterfaceParameters.SearchRadius, 
                                                rBinSearchInterfaceParameters.Results.begin(), 
                                                rBinSearchInterfaceParameters.ListOfDistances.begin(), 
                                                rBinSearchInterfaceParameters.NumberOfResults);                                     
    double minimum_distance = 1e14;

    KRATOS_ERROR_IF(obtained_results_2 == 0) << ":::[SnakeCutSbmProcess]::: No condition found in the skin model part for the surrogate node "
                                    << pSurrogateNode1->Id() << " at coordinates " << pSurrogateNode1->Coordinates() << std::endl;
    
    IndexType nearest_node_id = -1;
    for (IndexType k = 0; k < obtained_results_2; k++) {
        double current_distance = rBinSearchInterfaceParameters.ListOfDistances[k];   
        if (current_distance < minimum_distance) { 

            IndexType current_id = rBinSearchInterfaceParameters.Results[k]->Id();
            auto current_layers = rSkinSubModelPart.GetNode(current_id).GetValue(CONNECTED_LAYERS);

            bool layer_found = false;
            for (const auto& forced_layer : forced_layers) {
                if (std::find(current_layers.begin(), current_layers.end(), forced_layer) != current_layers.end()) {
                    layer_found = true;
                    break;
                }
            }
            if (!layer_found) continue;

            minimum_distance = current_distance;
            nearest_node_id = current_id;
        }
    }
    KRATOS_ERROR_IF(nearest_node_id == -1) << "::[SnakeSbmProcess]:: No suitable second condition found in search for projection of point: " <<
        p_surrogate_vertex_2 << " with first condition id: " << skin_node_id_1 << std::endl;

    KRATOS_ERROR_IF(obtained_results_2 == 0) << "::[SnakeSbmProcess]:: Zero points found in serch for projection of point: " <<
        p_surrogate_vertex_2 << std::endl;
    
    const IndexType id_closest_true_node_2 = nearest_node_id;

    // check if the first vertex was actually closer to this projection
    auto& r_projection_node = rSkinSubModelPart.GetNode(id_closest_true_node_2);

    if (norm_2(r_projection_node.Coordinates()-pSurrogateNode1->Coordinates()) < norm_2(r_projection_node.Coordinates()-pSurrogateNode2->Coordinates()))
    {
        pSurrogateNode1->SetValue(PROJECTION_NODE_ID, id_closest_true_node_2);

        auto connected_layers = rSkinSubModelPart.GetNode(id_closest_true_node_2).GetValue(CONNECTED_LAYERS);
        auto connected_conditions = rSkinSubModelPart.GetNode(id_closest_true_node_2).GetValue(CONNECTED_CONDITIONS);

        pSurrogateNode1->SetValue(CONNECTED_LAYERS, connected_layers);
        pSurrogateNode1->SetValue(CONNECTED_CONDITIONS, connected_conditions);
    }
    else
    {
        pSurrogateNode2->SetValue(PROJECTION_NODE_ID, id_closest_true_node_2);

        auto connected_layers = rSkinSubModelPart.GetNode(id_closest_true_node_2).GetValue(CONNECTED_LAYERS);
        auto connected_conditions = rSkinSubModelPart.GetNode(id_closest_true_node_2).GetValue(CONNECTED_CONDITIONS);

        pSurrogateNode2->SetValue(CONNECTED_LAYERS, connected_layers);
        pSurrogateNode2->SetValue(CONNECTED_CONDITIONS, connected_conditions);
    }
    return;
}

void SnakeCutSbmProcess::FindClosestTruePointToSurrogateVertexByNurbs()
{
    // // get the data
    // mEchoLevel = mThisParameters["echo_level"].GetInt();
    // std::string iga_model_part_name = mThisParameters["model_part_name"].GetString();
    // std::string skin_model_part_name = mThisParameters["skin_model_part_name"].GetString();
    // std::string skin_model_part_inner_initial_name = mThisParameters["skin_model_part_inner_initial_name"].GetString();
    
    // // Loop over the surrogate vertex and call the search in radius or something like that
    // // TODO: Outer boyndary
    // const bool is_inner = true;
    
    // std::string rSurrogateSubModelPartName; 
    // std::string skin_sub_model_part_name; 
    // rSurrogateSubModelPartName = "surrogate_inner";
    // skin_sub_model_part_name = "inner";

    // ModelPart& rCutSbmSubModelPart = mpIgaModelPart->CreateSubModelPart("extended_sbm");

    // ModelPart& rSkinSubModelPart = mpSkinModelPart->GetSubModelPart(skin_sub_model_part_name);
    // ModelPart& rSurrogateSubModelPart = mpIgaModelPart->GetSubModelPart(rSurrogateSubModelPartName);

    // auto p_surface = mpIgaModelPart->pGetGeometry(1);
    // IndexType id_brep_curve_on_surface = (mpIgaModelPart->GeometriesEnd()-1)->Id() + 1;
    
    // typedef NurbsSurfaceGeometry<3, PointerVector<NodeType>> NurbsSurfaceType;
    // auto p_nurbs_surface = std::dynamic_pointer_cast<NurbsSurfaceType>(
    //                         p_surface->pGetGeometryPart(Geometry<typename PointerVector<NodeType>::value_type>::BACKGROUND_GEOMETRY_INDEX));
    // IntegrationInfo surface_integration_info = p_nurbs_surface->GetDefaultIntegrationInfo();

    // // Get the mesh sizes from the surrogate model part
    // const Vector& knot_span_sizes = rSurrogateSubModelPart.GetParentModelPart().GetValue(KNOT_SPAN_SIZES);

    // double knot_span_reference_size = knot_span_sizes[0];
    // if (knot_span_sizes[1] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[1];}
    // if (knot_span_sizes.size() > 2) {if (knot_span_sizes[2] > knot_span_reference_size) {knot_span_reference_size = knot_span_sizes[2];}}

    // std::string condition_name = "ExtendedSbmSolidCondition";
    // std::string element_name = "ExtendedSbmSolidElement";
    // // std::string element_name = "SolidElement";

    // // Loop over the nodes of the surrogate sub model part
    // IndexType iel = 1;
    // SizeType number_of_shape_functions_derivatives = 5;

    // IndexType first_condition_id = rSurrogateSubModelPart.pGetElement(iel)->GetGeometry()[0].Id();
    // IndexType last_condition_id = rSurrogateSubModelPart.pGetElement(iel)->GetGeometry()[1].Id();

    // SizeType size_surrogate_loop = last_condition_id - first_condition_id + 1;

    // for (SizeType j = 0; j < size_surrogate_loop; ++j) {
    //     auto p_brep_geometry = mpIgaModelPart->pGetGeometry(6 + j);
    //     auto p_brep_curve_on_surface_surrogate1_surrogate2 = std::dynamic_pointer_cast<BrepCurveOnSurfaceType>(p_brep_geometry);

    //     KRATOS_ERROR_IF(!p_brep_curve_on_surface_surrogate1_surrogate2) <<  ":::[SnakeCutSbmProcess]::: the geometry with id " << p_brep_curve_on_surface_surrogate1_surrogate2->Id() 
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
    //         else if (skin_middle_point_local_coords[0] > rSkinSubModelPart.Nodes().front().Id())
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
    //         r_cut_sbm_sub_model_part, condition_name, id, PropertiesPointerType(), knot_span_sizes, surrogate_brep_middle_geometry, identifier_name);

        
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
