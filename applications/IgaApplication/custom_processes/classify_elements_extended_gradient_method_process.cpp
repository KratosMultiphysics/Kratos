//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Juan I. Camarotti

// Project includes
#include "classify_elements_extended_gradient_method_process.h"
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/geometry.h"

namespace Kratos
{

    ClassifyElementsExtendedGradientMethodProcess::ClassifyElementsExtendedGradientMethodProcess(
        Model& rModel, Parameters ThisParameters) : mpModel(&rModel), mParameters(ThisParameters)
    {
        ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    }

    void ClassifyElementsExtendedGradientMethodProcess::Execute(){

        ModelPart& background_mesh_model_part = mpModel->GetModelPart(mParameters["background_model_part_name"].GetString());
        ModelPart& embedded_body_model_part = mpModel->GetModelPart(mParameters["embedded_body_model_part_name"].GetString());
        bool keep_external_domain = mParameters["keep_external_domain"].GetBool();

        // Create sub model parts containing the elements intersected by the boundary and the active elements (to be integrated)
        ModelPart& intersected_elements_sub_model_part = background_mesh_model_part.CreateSubModelPart("intersected_elements");
        ModelPart& active_elements_sub_model_part = background_mesh_model_part.CreateSubModelPart("active_elements");
        ModelPart& interpolation_sub_model_part = background_mesh_model_part.CreateSubModelPart("interpolation");
        ModelPart& compute_error_sub_model_part = background_mesh_model_part.CreateSubModelPart("compute_error");

        // Create a vector containing the points defining the embedded geometry 
        std::vector<array_1d<double, 3>> polygon_vertices;
        for (auto condition_it = embedded_body_model_part.ConditionsBegin(); condition_it != embedded_body_model_part.ConditionsEnd(); condition_it++){
            GeometryType condition_geometry = condition_it->GetGeometry();

            // Get the points defining the line geometry
            array_1d<double,3> line_p1 = condition_geometry.GetPoint(0).Coordinates();
            polygon_vertices.push_back(line_p1);
        }
        
        // Recover the nurbs surface
        auto brep_surface = background_mesh_model_part.GeometriesEnd();
        for (auto geometry_it = background_mesh_model_part.GeometriesBegin(); geometry_it != background_mesh_model_part.GeometriesEnd(); geometry_it ++){
            std::string geometry_type = std::string(geometry_it->Info());
            if (geometry_type == "Brep surface"){
                brep_surface = geometry_it;
            }
        }

        // Get a pointer to the base class of the nurbs surface 
        auto nurbs_surface = brep_surface->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);

        // Get a pointer to the nurbs surface (derived class)
        auto p_nurbs_surface = dynamic_pointer_cast<NurbsSurfaceGeometry<3, PointerVector<Node>>>(nurbs_surface);

        // Get the knot vectors in the u and v directions
        std::vector<double> knot_vector_u, knot_vector_v;
        p_nurbs_surface->SpansLocalSpace(knot_vector_u, 0);
        p_nurbs_surface->SpansLocalSpace(knot_vector_v, 0);
        KRATOS_WATCH(knot_vector_u)

        // Resize the vector with the number of quadrature points inside the background mesh
        if (mpEntityPointVectorKDTree.size() != background_mesh_model_part.NumberOfElements()) {
            mpEntityPointVectorKDTree.resize(background_mesh_model_part.NumberOfElements());
        }

        // now fill the points vector
        IndexPartition<IndexType>(background_mesh_model_part.NumberOfElements()).for_each([&](const IndexType Index) {
            mpEntityPointVectorKDTree[Index] = Kratos::make_shared<EntityPoint<EntityType>>(*(background_mesh_model_part.ElementsBegin() + Index), Index);
        });

        // Constructing a search tree to efficiently look for integration points close to a knot span
        KDTree::Pointer pSearchTree =  Kratos::make_shared<ClassifyElementsExtendedGradientMethodProcess::KDTree>(mpEntityPointVectorKDTree.begin(), mpEntityPointVectorKDTree.end(), 100);
        // IndexType maximum_number_neighbours = std::pow((p_nurbs_surface->PolynomialDegreeU() + 1) * (p_nurbs_surface->PolynomialDegreeV() + 1), 2);
        IndexType maximum_number_neighbours = 500;

        std::cout << "Starting intersected elements classification" << std::endl;
        // Check if an element is intersected by the embedded geometry
        /* Idea:
            For every element in the parameter space:
            1) define if any of the lines representing the embedded boundary intersect the element
            2) define which quadrature points are inside the element 
        
        */
        for (IndexType i = 0; i < knot_vector_u.size() - 1; i++){
            for (IndexType j = 0; j < knot_vector_v.size() - 1; j++){
                // Element coordinates
                double u_min = knot_vector_u[i], u_max = knot_vector_u[i+1];
                double v_min = knot_vector_v[j], v_max = knot_vector_v[j+1];

                CoordinatesArrayType knot_span_center{0.5*(u_min + u_max), 0.5*(v_min + v_max), 0.0};
                double knot_span_diagonal = std::sqrt(std::pow(u_min - u_max, 2) + std::pow(v_min - v_max, 2));
                double radius = 3.0 * knot_span_diagonal;

                // Many lines could be inside a knot span. Thie vector saves the normal of each of this segments and then calculates an average normal
                std::vector<Vector> vector_unit_normal;
                
                IndexType id = 0;
                
                // Iterate over the conditions (each line is a condition)
                for (auto condition_it = embedded_body_model_part.ConditionsBegin(); condition_it != embedded_body_model_part.ConditionsEnd(); condition_it++){
                    GeometryType condition_geometry = condition_it->GetGeometry();

                    // Get the points defining the line geometry
                    array_1d<double,3> line_p1 = condition_geometry.GetPoint(0).Coordinates();
                    array_1d<double,3> line_p2 = condition_geometry.GetPoint(1).Coordinates();

                    // Verify intersection
                    bool intersects = LineIntersectsRectangle(line_p1, line_p2, u_min, u_max, v_min, v_max);
                    if (intersects == true){
                        Vector unit_normal = ComputeUnitNormal(line_p1, line_p2);
                        vector_unit_normal.push_back(unit_normal);

                        // If a knot span is intersected, look for the nearest integration point
                        NodeType::Pointer qp = Kratos::make_intrusive<NodeType>(
                            id, knot_span_center[0], knot_span_center[1], knot_span_center[2]
                        );
                        PointsArrayType knot_span_center_array;
                        knot_span_center_array.push_back(qp);
                        auto p_gp_geometry = Kratos::make_shared<Point3D<NodeType>>(
                            id, knot_span_center_array
                        );
                        LaplacianIGAElement knot_span_center_element{id, p_gp_geometry};
                        EntityPoint<EntityType> knot_span_center_entity_point(knot_span_center_element, id);
                        
                        // Create 2 vectors: one for the EntityPoint and other for the distances
                        std::vector<EntityPoint<EntityType>::Pointer> neighbours_entity_points;
                        std::vector<double> neighbours_entity_points_distances;
                        neighbours_entity_points.resize(maximum_number_neighbours);
                        neighbours_entity_points_distances.resize(maximum_number_neighbours);

                        // Looking for the m closest neighbours to the knot span center 
                        const auto number_of_neighbors = pSearchTree->SearchInRadius(
                                                            knot_span_center_entity_point,
                                                            radius,
                                                            neighbours_entity_points.begin(),
                                                            neighbours_entity_points_distances.begin(),
                                                            maximum_number_neighbours);

                        id += 1;
                        
                        for (IndexType k = 0; k < neighbours_entity_points.size(); k++){
                             if (neighbours_entity_points[k] != nullptr) { // Check if the pointer is not null
                                auto element = neighbours_entity_points[k]->GetEntity();
                                array_1d<double,3> quadrature_point_position = element.GetGeometry().Center();
                                IndexType quadrature_point_id = element.Id();
                                
                                if (IsPointInsideRectangle(quadrature_point_position[0], quadrature_point_position[1], u_min, u_max, v_min, v_max) == true){
                                    element.Set(ACTIVE, false);
                                    intersected_elements_sub_model_part.AddElement(background_mesh_model_part.pGetElement(quadrature_point_id));
                                }
                            }
                        }

                        // for (auto element_it = background_mesh_model_part.ElementsBegin(); element_it != background_mesh_model_part.ElementsEnd(); element_it ++){
                        //     array_1d<double,3> quadrature_point_position = element_it->GetGeometry().Center();

                        //     if (IsPointInsideRectangle(quadrature_point_position[0], quadrature_point_position[1], u_min, u_max, v_min, v_max) == true){
                        //         element_it->Set(ACTIVE, false);
                        //         intersected_elements_sub_model_part.AddElement(*(element_it.base()));
                        //         //element_it->SetValue(UNIT_NORMAL, unit_normal);
                        //     }
                        // }
                    }
                }

                Vector mean_unit_normal_vector = ZeroVector(3);
                for (auto normal_vector_it = vector_unit_normal.begin(); normal_vector_it != vector_unit_normal.end(); ++normal_vector_it) {
                    mean_unit_normal_vector[0] += (*normal_vector_it)[0];
                    mean_unit_normal_vector[1] += (*normal_vector_it)[1];
                    mean_unit_normal_vector[2] += (*normal_vector_it)[2];
                }

                // Divide each component by the number of vectors to compute the average
                size_t num_vectors = vector_unit_normal.size();
                if (num_vectors > 0) {
                    mean_unit_normal_vector[0] /= num_vectors;
                    mean_unit_normal_vector[1] /= num_vectors;
                    mean_unit_normal_vector[2] /= num_vectors;
                }
                // Normalize the normal vector 
                mean_unit_normal_vector /= norm_2(mean_unit_normal_vector);

                // Once the average normal inside a knot span is calculated, assign this value to the GP inside the knot span
                NodeType::Pointer qp = Kratos::make_intrusive<NodeType>(
                    id, knot_span_center[0], knot_span_center[1], knot_span_center[2]
                );
                PointsArrayType knot_span_center_array;
                knot_span_center_array.push_back(qp);
                auto p_gp_geometry = Kratos::make_shared<Point3D<NodeType>>(
                    id, knot_span_center_array
                );
                LaplacianIGAElement knot_span_center_element{id, p_gp_geometry};
                EntityPoint<EntityType> knot_span_center_entity_point(knot_span_center_element, id);

                // Create 2 vectors: one for the EntityPoint and other for the distances
                std::vector<EntityPoint<EntityType>::Pointer> neighbours_entity_points;
                std::vector<double> neighbours_entity_points_distances;
                neighbours_entity_points.resize(maximum_number_neighbours);
                neighbours_entity_points_distances.resize(maximum_number_neighbours);

                const auto number_of_neighbors = pSearchTree->SearchInRadius(
                    knot_span_center_entity_point,
                    radius,
                    neighbours_entity_points.begin(),
                    neighbours_entity_points_distances.begin(),
                    maximum_number_neighbours);         

                for (IndexType i = 0; i < neighbours_entity_points.size(); i++){
                    if (neighbours_entity_points[i] != nullptr) { // Check if the pointer is not null
                        auto element = neighbours_entity_points[i]->GetEntity();
                        array_1d<double,3> quadrature_point_position = element.GetGeometry().Center();
                        IndexType quadrature_point_id = element.Id();
                        
                        // If the GP is inside the knot span, assign the value of the normal vector
                        if (IsPointInsideRectangle(quadrature_point_position[0], quadrature_point_position[1], u_min, u_max, v_min, v_max) == true){
                            element.SetValue(UNIT_NORMAL, mean_unit_normal_vector);
                        }
                    }
                }   
            }
        }
        std::cout << "Finishing intersected elements classification" << std::endl;

        std::cout << "Starting inside elements classification" << std::endl;
        // Check if an element is inside the geometry using the ray tracing algorithm 
        for (IndexType i = 0; i < knot_vector_u.size() - 1; i++){
            for (IndexType j = 0; j < knot_vector_v.size() - 1; j++){
                // Element coordinates
                double u_min = knot_vector_u[i], u_max = knot_vector_u[i+1];
                double v_min = knot_vector_v[j], v_max = knot_vector_v[j+1];

                // Set the center and radius of the kd tree search 
                CoordinatesArrayType knot_span_center{0.5*(u_min + u_max), 0.5*(v_min + v_max), 0.0};
                double knot_span_diagonal = std::sqrt(std::pow(u_min - u_max, 2) + std::pow(v_min - v_max, 2));
                double radius = 1.3 * knot_span_diagonal;

                // Define the four corners of the rectangle
                array_1d<double, 3> p1{u_min, v_min, 0.0};  // Corner 1
                array_1d<double, 3> p2{u_max, v_min, 0.0};  // Corner 2
                array_1d<double, 3> p3{u_max, v_max, 0.0};  // Corner 3
                array_1d<double, 3> p4{u_min, v_max, 0.0};  // Corner 4

                // Check if the rectangle is inside or outside the polygon using the ray tracing algorithm 
                auto is_rectangle_inside = IsRectangleInsidePolygon(p1, p2, p3, p4, polygon_vertices);
                auto number_of_points_inside = is_rectangle_inside.second;

                // Creating an entity point in the knot span middle 
                NodeType::Pointer qp = Kratos::make_intrusive<NodeType>(
                    0, knot_span_center[0], knot_span_center[1], knot_span_center[2]
                );
                PointsArrayType knot_span_center_array;
                knot_span_center_array.push_back(qp);
                auto p_gp_geometry = Kratos::make_shared<Point3D<NodeType>>(
                    0, knot_span_center_array
                );
                LaplacianIGAElement knot_span_center_element{0, p_gp_geometry};
                EntityPoint<EntityType> knot_span_center_entity_point(knot_span_center_element, 0);

                // Create 2 vectors: one for the EntityPoint and other for the distances
                std::vector<EntityPoint<EntityType>::Pointer> neighbours_entity_points;
                std::vector<double> neighbours_entity_points_distances;
                neighbours_entity_points.resize(maximum_number_neighbours);
                neighbours_entity_points_distances.resize(maximum_number_neighbours);

                // Perform the search in radius for looking the m closest neighbours
                const auto number_of_neighbors = pSearchTree->SearchInRadius(
                                                    knot_span_center_entity_point,
                                                    radius,
                                                    neighbours_entity_points.begin(),
                                                    neighbours_entity_points_distances.begin(),
                                                    maximum_number_neighbours);
                
                // If we want to keep the knot spans outside the embedded geometry, deactivate those elements which are inside
                if (is_rectangle_inside.first == true && keep_external_domain == true){
                    for (IndexType k = 0; k < neighbours_entity_points.size(); k++){
                        if (neighbours_entity_points[k] != nullptr) { // Check if the pointer is not null
                            auto element = neighbours_entity_points[k]->GetEntity();
                            array_1d<double,3> quadrature_point_position = element.GetGeometry().Center();
                            IndexType quadrature_point_id = element.Id();
                                    
                            if (IsPointInsideRectangle(quadrature_point_position[0], quadrature_point_position[1], u_min, u_max, v_min, v_max) == true){
                                background_mesh_model_part.pGetElement(quadrature_point_id)->Set(ACTIVE, false);
                            }
                        }
                    }   
                    // for (auto element_it = background_mesh_model_part.ElementsBegin(); element_it != background_mesh_model_part.ElementsEnd(); element_it ++){
                    //     array_1d<double,3> quadrature_point_position = element_it->GetGeometry().Center();
                        
                    //     if (IsPointInsideRectangle(quadrature_point_position[0], quadrature_point_position[1], u_min, u_max, v_min, v_max) == true){
                    //         element_it->Set(ACTIVE, false);
                    //     }
                    // }
                }
                else if (is_rectangle_inside.first == false && keep_external_domain == true){
                    for (IndexType k = 0; k < neighbours_entity_points.size(); k++){
                        if (neighbours_entity_points[k] != nullptr) { // Check if the pointer is not null
                            auto element = neighbours_entity_points[k]->GetEntity();
                            array_1d<double,3> quadrature_point_position = element.GetGeometry().Center();
                            IndexType quadrature_point_id = element.Id();
                            
                            // If the rectangle is outside and we keep the external domain, we add to the corresponding GP to the active elements
                            if (IsPointInsideRectangle(quadrature_point_position[0], quadrature_point_position[1], u_min, u_max, v_min, v_max) == true){
                                active_elements_sub_model_part.AddElement(background_mesh_model_part.pGetElement(quadrature_point_id));
                            }

                            if (IsPointInsideRectangle(quadrature_point_position[0], quadrature_point_position[1], u_min, u_max, v_min, v_max) == true && number_of_points_inside == 0){
                                interpolation_sub_model_part.AddElement(background_mesh_model_part.pGetElement(quadrature_point_id));
                            }

                            if (IsPointInsideRectangle(quadrature_point_position[0], quadrature_point_position[1], u_min, u_max, v_min, v_max) == true && number_of_points_inside == 0){
                                for (auto element_it_2 = intersected_elements_sub_model_part.ElementsBegin(); element_it_2 != intersected_elements_sub_model_part.ElementsEnd(); element_it_2++){
                                    array_1d<double,3> gp_position_intersected_elements = element_it_2->GetGeometry().Center();

                                    double distance = norm_2(gp_position_intersected_elements - quadrature_point_position);

                                    if (distance > 0.05){
                                        compute_error_sub_model_part.AddElement(background_mesh_model_part.pGetElement(quadrature_point_id));
                                    }
                                }
                            }
                        }
                    }   
                    // for (auto element_it = background_mesh_model_part.ElementsBegin(); element_it != background_mesh_model_part.ElementsEnd(); element_it ++){
                    //     array_1d<double,3> quadrature_point_position = element_it->GetGeometry().Center();

                    //         if (IsPointInsideRectangle(quadrature_point_position[0], quadrature_point_position[1], u_min, u_max, v_min, v_max) == true){
                    //             active_elements_sub_model_part.AddElement(*(element_it.base()));
                    //         }

                    //         if (IsPointInsideRectangle(quadrature_point_position[0], quadrature_point_position[1], u_min, u_max, v_min, v_max) == true && number_of_points_inside == 0){
                    //             interpolation_sub_model_part.AddElement(*(element_it.base()));
                    //         }

                    //         if (IsPointInsideRectangle(quadrature_point_position[0], quadrature_point_position[1], u_min, u_max, v_min, v_max) == true && number_of_points_inside == 0){
                    //             for (auto element_it_2 = intersected_elements_sub_model_part.ElementsBegin(); element_it_2 != intersected_elements_sub_model_part.ElementsEnd(); element_it_2++){
                    //                 array_1d<double,3> gp_position_intersected_elements = element_it_2->GetGeometry().Center();

                    //                 double distance = norm_2(gp_position_intersected_elements - quadrature_point_position);

                    //                 if (distance > 0.0){
                    //                     compute_error_sub_model_part.AddElement(*(element_it_2.base()));
                    //                 }
                    //             }
                    //         }
                        
                    // }
                }
            }
        }
        std::cout << "Finishing inside elements classification" << std::endl;

        // For every segment defining the embedded geometry, create a brep curve on surface geometry and create conditions with integration points
        for (auto condition_it = embedded_body_model_part.ConditionsBegin(); condition_it != embedded_body_model_part.ConditionsEnd(); condition_it++){
            GeometryType condition_geometry = condition_it->GetGeometry();
            
            Geometry<GeometricalObject::NodeType>::PointsArrayType brep_curve_points;
            brep_curve_points.push_back(embedded_body_model_part.pGetNode(condition_geometry.GetPoint(0).Id()));
            brep_curve_points.push_back(embedded_body_model_part.pGetNode(condition_geometry.GetPoint(1).Id()));
                                
            // Creating the brep curve knot vector (simplest knot vector for a line segment)
            Vector knot_vector = ZeroVector(4);
            knot_vector[0] = 0.0;
            knot_vector[1] = 0.0;
            knot_vector[2] = 1.0;
            knot_vector[3] = 1.0;
                                
            // Setting the polinomial degree of the b-rep curve 
            int p = 1;

            // Creating a pointer to the b-rep curve and to the nurbs curve on surface
            NurbsCurveGeometry<2, PointerVector<NodeType>>::Pointer p_brep_curve = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<NodeType>>>(brep_curve_points, p, knot_vector);
            NurbsCurveOnSurfaceGeometry<3, PointerVector<NodeType>, PointerVector<NodeType>> nurbs_curve_on_surface(p_nurbs_surface, p_brep_curve);

            // Create integration points 
            IntegrationPointsArrayType integration_points;
            IntegrationPointsArrayType integration_points_modified;
            IntegrationInfo integration_info = nurbs_curve_on_surface.GetDefaultIntegrationInfo();
            nurbs_curve_on_surface.CreateIntegrationPoints(integration_points, integration_info);
            
            // Do we really need this? 
            for (auto int_point_iterator = integration_points.begin(); int_point_iterator != integration_points.end(); int_point_iterator++){
                if (int_point_iterator->Weight() > 5e-3){
                    integration_points_modified.push_back(*int_point_iterator);
                }
            }
            
            // Create quadrature point geometries
            GeometriesArrayType quadrature_points_geometries(integration_points.size()); 
            IndexType number_of_shape_functions_derivatives = 2;
            nurbs_curve_on_surface.CreateQuadraturePointGeometries(quadrature_points_geometries, number_of_shape_functions_derivatives, integration_points_modified, integration_info);
            
            // For each quadrature point, create a new condition
            IndexType size = quadrature_points_geometries.size();
            IndexType i = 0;
            const IndexType last_condition_id = (intersected_elements_sub_model_part.GetRootModelPart().NumberOfConditions() == 0)
                                         ? 1
                                         : (intersected_elements_sub_model_part.GetRootModelPart().ConditionsEnd() - 1)->Id() + 1 ;
            for (auto element_it = quadrature_points_geometries.begin(); element_it != quadrature_points_geometries.end(); element_it++){
                intersected_elements_sub_model_part.AddCondition(Kratos::make_intrusive<Condition>(
                        last_condition_id + i, quadrature_points_geometries(i)));
                i += 1;
            }
        }

        for (auto element_it = compute_error_sub_model_part.ElementsBegin(); element_it != compute_error_sub_model_part.ElementsEnd(); element_it++){
            auto gauss_point_weight = element_it->GetGeometry().IntegrationPoints()[0].Weight();
            Vector weights(1); 
            weights[0] = gauss_point_weight;
            element_it->SetValue(INTEGRATION_WEIGHTS, weights);
        }
    // KRATOS_WATCH("outside classification")
    }

    Vector ClassifyElementsExtendedGradientMethodProcess::ComputeUnitNormal(const array_1d<double,3>& line_p1, const array_1d<double,3>& line_p2){
        Vector normal(3);
        normal[0] = - (line_p2[1] - line_p1[1]);
        normal[1] = (line_p2[0] - line_p1[0]);
        normal[2] = 0.0;
        return - normal/norm_2(normal);
    }

    // This method checks if a point is inside a polygon using the ray tracing algorithm
    bool ClassifyElementsExtendedGradientMethodProcess::IsPointInsidePolygon(const array_1d<double, 3>& point, const std::vector<array_1d<double, 3>>& polygon_vertices) {
        int n = polygon_vertices.size();
        bool inside = false;

        // Loop over each edge of the polygon
        for (int i = 0, j = n - 1; i < n; j = i++) {
            const array_1d<double, 3>& p1 = polygon_vertices[i];
            const array_1d<double, 3>& p2 = polygon_vertices[j];

            // Check if the ray from point crosses the edge from p1 to p2
            if (((p1[1] > point[1]) != (p2[1] > point[1])) && 
                (point[0] < (p2[0] - p1[0]) * (point[1] - p1[1]) / (p2[1] - p1[1]) + p1[0])) {
                inside = !inside;  // Toggle the inside flag
            }
        }

        return inside;  // Return whether the point is inside the polygon
    }

    std::pair<bool, int>  ClassifyElementsExtendedGradientMethodProcess::IsRectangleInsidePolygon(const array_1d<double, 3>& p1, const array_1d<double, 3>& p2,
                               const array_1d<double, 3>& p3, const array_1d<double, 3>& p4, std::vector<array_1d<double, 3>> polygon_vertices) {
        
        bool is_p1_inside = IsPointInsidePolygon(p1, polygon_vertices);
        bool is_p2_inside = IsPointInsidePolygon(p2, polygon_vertices);
        bool is_p3_inside = IsPointInsidePolygon(p3, polygon_vertices);
        bool is_p4_inside = IsPointInsidePolygon(p4, polygon_vertices);

        int number_of_points_inside = is_p1_inside + is_p2_inside + is_p3_inside + is_p4_inside;

        // Check if all the rectangle corners are inside the polygon
        return std::make_pair(is_p1_inside &&
           is_p2_inside &&
           is_p3_inside &&
           is_p4_inside, number_of_points_inside);
    }

    // Function to verify if a point is inside the rectangle
    bool ClassifyElementsExtendedGradientMethodProcess::IsPointInsideRectangle(double x, double y, double u_min, double u_max, double v_min, double v_max) {
        return (x > u_min && x < u_max && y > v_min && y < v_max);
    }

    // Function to verify if 2 segments intersect
    bool ClassifyElementsExtendedGradientMethodProcess::DoSegmentsIntersect(
        const array_1d<double,3>& p1, const array_1d<double,3>& p2,
        const array_1d<double,3>& q1, const array_1d<double,3>& q2) 
    {
        
        auto Orientation = [](const array_1d<double,3>& a, const array_1d<double,3>& b, const array_1d<double,3>& c) {
            return (b[1] - a[1]) * (c[0] - b[0]) - (b[0] - a[0]) * (c[1] - b[1]);
        };

        double o1 = Orientation(p1, p2, q1);
        double o2 = Orientation(p1, p2, q2);
        double o3 = Orientation(q1, q2, p1);
        double o4 = Orientation(q1, q2, p2);

        // Comprobar si las orientaciones se cruzan
        if (o1 * o2 < 0 && o3 * o4 < 0) return true;

        return false; // No se intersectan
    }

    // Main function to verify the line-rectangle intersection
    bool ClassifyElementsExtendedGradientMethodProcess::LineIntersectsRectangle(
        const array_1d<double,3>& line_p1, 
        const array_1d<double,3>& line_p2,
        double u_min, double u_max, double v_min, double v_max) 
    {
        // Verify if any of the points defining the line is inside the rectangle
        if (IsPointInsideRectangle(line_p1[0], line_p1[1], u_min, u_max, v_min, v_max) ||
            IsPointInsideRectangle(line_p2[0], line_p2[1], u_min, u_max, v_min, v_max)) 
        {
            return true;
        }

        // Define the boundary of the rectangle as segments
        array_1d<double, 3> p1{u_min, v_min, 0.0};
        array_1d<double, 3> p2{u_max, v_min, 0.0};
        array_1d<double, 3> p3{u_max, v_max, 0.0};
        array_1d<double, 3> p4{u_min, v_max, 0.0};

        // Using a vector of pairs of points to represent the edges
        std::vector<std::pair<array_1d<double, 3>, array_1d<double, 3>>> rect_edges = {
            {p1, p2},  // Segment 1
            {p2, p3},  // Segment 2
            {p3, p4},  // Segment 3
            {p4, p1}   // Segment 4
        };

        // Verify intersection with each side of the rectangle
        for (const auto& edge : rect_edges) {
            if (DoSegmentsIntersect(line_p1, line_p2, edge.first, edge.second)) {
                return true;
            }
        }

        return false; // No hay intersección
    }

} // End namespace Kratos
