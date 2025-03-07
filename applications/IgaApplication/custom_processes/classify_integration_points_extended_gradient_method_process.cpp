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
#include "classify_integration_points_extended_gradient_method_process.h"
#include "geometries/nurbs_volume_geometry.h"
#include "geometries/geometry.h"

namespace Kratos
{

    ClassifyIntegrationPointsExtendedGradientMethodProcess::ClassifyIntegrationPointsExtendedGradientMethodProcess(
        Model& rModel, Parameters ThisParameters) : mpModel(&rModel), mParameters(ThisParameters)
    {
        ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    }


    void ClassifyIntegrationPointsExtendedGradientMethodProcess::Execute() {
        mpBackgroundMeshModelPart = &mpModel->GetModelPart(mParameters["background_model_part_name"].GetString());
        mpSkinOuterLoopModelPart = &mpModel->GetModelPart(mParameters["skin_model_part_outer_name"].GetString());
        mNumberOfInnerLoops = mParameters["number_of_inner_loops"].GetInt();
        if (mNumberOfInnerLoops > 0) {
            mpSkinInnerLoopModelPart = &mpModel->GetModelPart(mParameters["skin_model_part_inner_name"].GetString());
        }

        // Create sub model parts containing the elements intersected by the boundary and the active elements
        ModelPart& intersected_elements_sub_model_part = mpBackgroundMeshModelPart->CreateSubModelPart("intersected_elements");
        ModelPart& active_elements_sub_model_part = mpBackgroundMeshModelPart->CreateSubModelPart("active_elements");

        // Create a vector containing the points defining the outer loop
        for (auto condition_it_outer = mpSkinOuterLoopModelPart->ConditionsBegin(); condition_it_outer != mpSkinOuterLoopModelPart->ConditionsEnd(); condition_it_outer++){
            GeometryType condition_geometry = condition_it_outer->GetGeometry();

            // Get the points defining the line geometry
            array_1d<double,3> line_p1 = condition_geometry.GetPoint(0).Coordinates();
            mOuterLoopPolygon.push_back(line_p1);
        }

        // Create a vector containing the points defining the inner loop
        if (mNumberOfInnerLoops > 0){
            for (auto condition_it_inner = mpSkinInnerLoopModelPart->ConditionsBegin(); condition_it_inner != mpSkinInnerLoopModelPart->ConditionsEnd(); condition_it_inner++){
                GeometryType condition_geometry = condition_it_inner->GetGeometry();

                // Get the points defining the line geometry
                array_1d<double,3> line_p1 = condition_geometry.GetPoint(0).Coordinates();
                mInnerLoopPolygon.push_back(line_p1);
            }
        }

        // Store elements to be added later (avoiding iterator invalidation)
        std::vector<Element::Pointer> elements_to_intersect;
        std::vector<Element::Pointer> elements_to_activate;

        // Recover the background NURBS surface
        auto brep_surface = mpBackgroundMeshModelPart->GeometriesEnd();
        for (auto geometry_it = mpBackgroundMeshModelPart->GeometriesBegin(); geometry_it != mpBackgroundMeshModelPart->GeometriesEnd(); ++geometry_it) {
            std::string geometry_type = std::string(geometry_it->Info());
            if (geometry_type == "Brep surface") {
                brep_surface = geometry_it;
            }
        }

        auto nurbs_surface = brep_surface->pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX);
        auto p_nurbs_surface = dynamic_pointer_cast<NurbsSurfaceGeometry<3, PointerVector<Node>>>(nurbs_surface);

        std::vector<double> knot_vector_u, knot_vector_v;
        p_nurbs_surface->SpansLocalSpace(knot_vector_u, 0);
        p_nurbs_surface->SpansLocalSpace(knot_vector_v, 0);

        double polynomial_degree_u = p_nurbs_surface->PolynomialDegreeU();
        double polynomial_degree_v = p_nurbs_surface->PolynomialDegreeV();

        // Classification procedure
        for (auto integration_point_it = mpBackgroundMeshModelPart->ElementsBegin(); integration_point_it != mpBackgroundMeshModelPart->ElementsEnd();) {
            auto integration_point_center = integration_point_it->pGetGeometry()->Center();

            bool found_knot_span = false;
            double u_min = 0.0, u_max = 0.0;
            double v_min = 0.0, v_max = 0.0;

            for (IndexType i = 0; i < knot_vector_u.size() - 1; i++) {
                for (IndexType j = 0; j < knot_vector_v.size() - 1; j++) {
                    if (integration_point_center.X() > knot_vector_u[i] && integration_point_center.X() < knot_vector_u[i + 1] &&
                        integration_point_center.Y() > knot_vector_v[j] && integration_point_center.Y() < knot_vector_v[j + 1]) {
                        u_min = knot_vector_u[i];
                        u_max = knot_vector_u[i + 1];
                        v_min = knot_vector_v[j];
                        v_max = knot_vector_v[j + 1];

                        found_knot_span = true;
                        break;
                    }
                }
                if (found_knot_span) {
                    break;
                }
            }

            Vector knot_span_normal_vector = ZeroVector(3);
            bool is_knot_span_intersected = IsKnotSpanIntersected(u_min, u_max, v_min, v_max, knot_span_normal_vector);

            size_t num_integration_points = (polynomial_degree_u + 1) * (polynomial_degree_v + 1);

            // Advance the end iterator safely
            auto integration_point_it_end = integration_point_it;
            for (size_t i = 0; i < num_integration_points && integration_point_it_end != mpBackgroundMeshModelPart->ElementsEnd(); ++i) {
                ++integration_point_it_end;
            }

            if (is_knot_span_intersected) {
                while (integration_point_it != integration_point_it_end) {
                    integration_point_it->Set(ACTIVE, false);
                    elements_to_intersect.push_back(integration_point_it->shared_from_this());
                    integration_point_it->SetValue(UNIT_NORMAL, knot_span_normal_vector);
                    ++integration_point_it;
                }
                continue;
            }

            bool is_knot_span_active = IsKnotSpanActive(u_min, u_max, v_min, v_max);

            if (is_knot_span_active) {
                while (integration_point_it != integration_point_it_end) {
                    integration_point_it->Set(ACTIVE, true);
                    elements_to_activate.push_back(integration_point_it->shared_from_this());
                    ++integration_point_it;
                }
                continue;
            }

            if (!is_knot_span_active) {
                while (integration_point_it != integration_point_it_end) {
                    integration_point_it->Set(ACTIVE, false);
                    ++integration_point_it;
                }
                continue;
            }
        }

        // Now modify the model safely after iteration
        for (auto& elem : elements_to_intersect) {
            intersected_elements_sub_model_part.AddElement(elem);
        }
        for (auto& elem : elements_to_activate) {
            active_elements_sub_model_part.AddElement(elem);
        }

        // Create the integration points in the outer and inner loops
        CreateBrepCurveOnSurfaceIntegrationPoints(p_nurbs_surface, intersected_elements_sub_model_part, mpSkinOuterLoopModelPart);
        if (mNumberOfInnerLoops > 0) {
            CreateBrepCurveOnSurfaceIntegrationPoints(p_nurbs_surface, intersected_elements_sub_model_part, mpSkinInnerLoopModelPart);
        }
    }

    bool ClassifyIntegrationPointsExtendedGradientMethodProcess::IsKnotSpanIntersected(
        const double u_min, const double u_max, const double v_min, const double v_max, Vector& unit_normal_vector){
        
        std::vector<Vector> vector_of_unit_normals, vector_of_unit_normals_inner, vector_of_unit_normals_outer;

        // Iterate over the lines in the outer loop and check if a knot span is intersected
        for (auto condition_it_outer_loop = mpSkinOuterLoopModelPart->ConditionsBegin(); condition_it_outer_loop != mpSkinOuterLoopModelPart->ConditionsEnd(); condition_it_outer_loop++){
            GeometryType condition_geometry = condition_it_outer_loop->GetGeometry();

            // Get the points defining the line geometry
            array_1d<double,3> line_p1 = condition_geometry.GetPoint(0).Coordinates();
            array_1d<double,3> line_p2 = condition_geometry.GetPoint(1).Coordinates();

            // Verify intersection
            bool intersects = LineIntersectsRectangle(line_p1, line_p2, u_min, u_max, v_min, v_max);
            
            if (intersects == true){
                Vector unit_normal = ComputeUnitNormal(line_p1, line_p2);
                vector_of_unit_normals_outer.push_back(unit_normal);
            }
        }

        // Iterate over the lines in the inner loop and check if a knot span is intersected
        if (mNumberOfInnerLoops > 0){
            for (auto condition_it_inner_loop = mpSkinInnerLoopModelPart->ConditionsBegin(); condition_it_inner_loop != mpSkinInnerLoopModelPart->ConditionsEnd(); condition_it_inner_loop++){
                GeometryType condition_geometry = condition_it_inner_loop->GetGeometry();

                // Get the points defining the line geometry
                array_1d<double,3> line_p1 = condition_geometry.GetPoint(0).Coordinates();
                array_1d<double,3> line_p2 = condition_geometry.GetPoint(1).Coordinates();

                // Verify intersection
                bool intersects = LineIntersectsRectangle(line_p1, line_p2, u_min, u_max, v_min, v_max);
                
                if (intersects == true){
                    Vector unit_normal = ComputeUnitNormal(line_p1, line_p2);
                    vector_of_unit_normals_inner.push_back(unit_normal);
                }
            }
        }

        KRATOS_ERROR_IF(vector_of_unit_normals_outer.size() > 0 && vector_of_unit_normals_inner.size() > 0) << "::[ClassifyIntegrationPointsExtendedGradientMethodProcess]::" 
                << "A knot span cannot be simultaneously intersected by an inner and outer loop. Decrease the knot span size." << std::endl;

        if (vector_of_unit_normals_outer.size() > 0) vector_of_unit_normals = vector_of_unit_normals_outer;
        if (vector_of_unit_normals_inner.size() > 0) vector_of_unit_normals = vector_of_unit_normals_inner;

        size_t number_of_segments_intersecting_knot_span = vector_of_unit_normals.size();

        // If the number of segments intersecting a knot span is zero, return that the knot span is not intersected
        if (number_of_segments_intersecting_knot_span == 0) return false;
        
        // Calculate an average unit normal vector for the intersected knot spans 
        Vector average_unit_normal_vector = ZeroVector(3);
        for (auto normal_vector_it = vector_of_unit_normals.begin(); normal_vector_it != vector_of_unit_normals.end(); ++normal_vector_it) {
            average_unit_normal_vector[0] += (*normal_vector_it)[0];
            average_unit_normal_vector[1] += (*normal_vector_it)[1];
            average_unit_normal_vector[2] += (*normal_vector_it)[2];
        }

        // Divide each component by the number of vectors to compute the average
        if (number_of_segments_intersecting_knot_span > 0) {
            average_unit_normal_vector[0] /= number_of_segments_intersecting_knot_span;
            average_unit_normal_vector[1] /= number_of_segments_intersecting_knot_span;
            average_unit_normal_vector[2] /= number_of_segments_intersecting_knot_span;
        }
        // Normalize the normal vector 
        average_unit_normal_vector /= norm_2(average_unit_normal_vector);
        
        // Return that the knot span is intersected ant the unit normal vector
        unit_normal_vector = average_unit_normal_vector;
        return true;
    }

    bool ClassifyIntegrationPointsExtendedGradientMethodProcess::IsKnotSpanActive(const double u_min, const double u_max, const double v_min, const double v_max){
        // For an element to be active, it needs to be inside the outer loop and outside the inner loop 
        // Define the four corners of the rectangle
        array_1d<double, 3> p1{u_min, v_min, 0.0}; 
        array_1d<double, 3> p2{u_max, v_min, 0.0};  
        array_1d<double, 3> p3{u_max, v_max, 0.0};  
        array_1d<double, 3> p4{u_min, v_max, 0.0};  

        // Check first if the knot span is inside or outside the outer polygon using the ray tracing algorithm
        auto is_rectangle_inside_outer_loop = IsRectangleInsidePolygon(p1, p2, p3, p4, mOuterLoopPolygon);

        //  If the knot span is outside the outer loop, return 0
        if (is_rectangle_inside_outer_loop.first == 0) return false; 
        if (is_rectangle_inside_outer_loop.first == 1 && mNumberOfInnerLoops == 0) return true; 
        
        if(mNumberOfInnerLoops > 0){
            // Check if the knot span is inside or outside the inner loop
            auto is_rectangle_inside_inner_loop = IsRectangleInsidePolygon(p1, p2, p3, p4, mInnerLoopPolygon);

            // If the knot span is inside the outer loop and inside the inner loop, return 0
            if (is_rectangle_inside_outer_loop.first == 1 && is_rectangle_inside_inner_loop.first == 1) return false; 

            // If the knot span is inside the outer loop and outside the inner loop, return 1
            if(is_rectangle_inside_outer_loop.first == 1 && is_rectangle_inside_inner_loop.first == 0) return true;
        }
    }

    // Main function to verify the line-rectangle intersection
    bool ClassifyIntegrationPointsExtendedGradientMethodProcess::LineIntersectsRectangle(
        const array_1d<double,3>& line_p1, 
        const array_1d<double,3>& line_p2,
        double u_min, double u_max, double v_min, double v_max) 
    {
        // Intersetion points between the line and the knot span
        std::vector<array_1d<double,3>> intersection_points;

        // Verify if any of the points defining the line is inside the rectangle
        bool is_p1_inside = IsPointInsideRectangle(line_p1[0], line_p1[1], u_min, u_max, v_min, v_max);
        bool is_p2_inside = IsPointInsideRectangle(line_p2[0], line_p2[1], u_min, u_max, v_min, v_max);

        if (is_p1_inside) {
            array_1d<double,3> intersection_point;
            intersection_point[0] = line_p1[0];
            intersection_point[1] = line_p1[1];
            intersection_point[2] = line_p1[2];
            intersection_points.push_back(intersection_point);
        }
        if (is_p2_inside) {
            array_1d<double,3> intersection_point;
            intersection_point[0] = line_p2[0];
            intersection_point[1] = line_p2[1];
            intersection_point[2] = line_p2[2];
            intersection_points.push_back(intersection_point);
        }

        // If no point is inside the rectangle, check if the line intersects any segment of the knot spans
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

        double epsilon = 1e-10;
        for (auto& edge : rect_edges) {
            array_1d<double,3> condition_knot_span_edge_intersection;
            if (DoSegmentsIntersect(line_p1, line_p2, edge.first, edge.second, condition_knot_span_edge_intersection)) {
                // Check if the intersection point already exists
                bool is_duplicate = false;
                for (const auto& existing_point : intersection_points) {
                    if (std::abs(existing_point[0] - condition_knot_span_edge_intersection[0]) < epsilon &&
                        std::abs(existing_point[1] - condition_knot_span_edge_intersection[1]) < epsilon) {
                        is_duplicate = true;
                        break;
                    }
                }

                // Only add unique intersection points
                if (!is_duplicate) {
                    intersection_points.push_back(condition_knot_span_edge_intersection);
                }
            }
        }
        
        // Calculate the distances from the intersection points to the knot span vertices
        if (intersection_points.size() == 2){
            array_1d<double, 3> distance_intersection_point_1_to_p1 = p1 - intersection_points[0];
            array_1d<double, 3> distance_intersection_point_1_to_p2 = p2 - intersection_points[0];
            array_1d<double, 3> distance_intersection_point_1_to_p3 = p3 - intersection_points[0];
            array_1d<double, 3> distance_intersection_point_1_to_p4 = p4 - intersection_points[0];
            array_1d<double, 3> distance_intersection_point_2_to_p1 = p1 - intersection_points[1];
            array_1d<double, 3> distance_intersection_point_2_to_p2 = p2 - intersection_points[1];
            array_1d<double, 3> distance_intersection_point_2_to_p3 = p3 - intersection_points[1];
            array_1d<double, 3> distance_intersection_point_2_to_p4 = p4 - intersection_points[1];

            // Check if the intersections lie on the knot span sides
            double u_product_int_point_1 = 
                (distance_intersection_point_1_to_p1[0] * distance_intersection_point_1_to_p2[0] * 
                distance_intersection_point_1_to_p3[0] * distance_intersection_point_1_to_p4[0]) / (norm_2(intersection_points[0]-intersection_points[1]));

            double u_product_int_point_2 = 
                 (distance_intersection_point_2_to_p1[0] * distance_intersection_point_2_to_p2[0] * 
                distance_intersection_point_2_to_p3[0] * distance_intersection_point_2_to_p4[0]) / (norm_2(intersection_points[0]-intersection_points[1]));

            double v_product_int_point_1 = 
                (distance_intersection_point_1_to_p1[1] * distance_intersection_point_1_to_p2[1] * 
                distance_intersection_point_1_to_p3[1] * distance_intersection_point_1_to_p4[1]) / (norm_2(intersection_points[0]-intersection_points[1]));

            double v_product_int_point_2 = 
                (distance_intersection_point_2_to_p1[1] * distance_intersection_point_2_to_p2[1] * 
                distance_intersection_point_2_to_p3[1] * distance_intersection_point_2_to_p4[1]) / (norm_2(intersection_points[0]-intersection_points[1]));

            if (u_product_int_point_1 == 0 && u_product_int_point_2 == 0 && v_product_int_point_1 == 0 && v_product_int_point_2 == 0) {
                return false;
            }
        }

        // Check if we have exactly two intersections at the same vertex, and ignore it
        if (intersection_points.size() == 2) {
            if (std::abs(intersection_points[0][0] - intersection_points[1][0]) < epsilon &&
                std::abs(intersection_points[0][1] - intersection_points[1][1]) < epsilon && 
                std::abs(intersection_points[0][2] - intersection_points[1][2]) < epsilon){
                return false;
            }
            return true;
        }

        return false; // There is no intersection
    }

     // Function to verify if a point is inside the rectangle
    bool ClassifyIntegrationPointsExtendedGradientMethodProcess::IsPointInsideRectangle(double x, double y, double u_min, double u_max, double v_min, double v_max){
        if (x >= u_min && x <= u_max && y >= v_min && y <= v_max){
            return true;
        }
        return false;
    }

     // Function to verify if 2 segments intersect
    bool ClassifyIntegrationPointsExtendedGradientMethodProcess::DoSegmentsIntersect(
        const array_1d<double,3>& p1, const array_1d<double,3>& p2,
        const array_1d<double,3>& q1, const array_1d<double,3>& q2,
        array_1d<double,3>& line_knot_span_edge_intersection_point) 
    {
        double epsilon = 1e-7; // Small tolerance for numerical precision

        auto Orientation = [epsilon](const array_1d<double,3>& a, const array_1d<double,3>& b, const array_1d<double,3>& c) {
            double val = (b[1] - a[1]) * (c[0] - b[0]) - (b[0] - a[0]) * (c[1] - b[1]);
            return (std::abs(val) < epsilon) ? 0 : (val > 0 ? 1 : -1);
        };

        // Lambda function to check if a point p is on a given segment defined by a and b
        auto IsPointOnSegment = [](const array_1d<double,3>& a, const array_1d<double,3>& b, const array_1d<double,3>& p) {
            return std::min(a[0], b[0]) <= p[0] && p[0] <= std::max(a[0], b[0]) &&
                std::min(a[1], b[1]) <= p[1] && p[1] <= std::max(a[1], b[1]);
        };

        double o1 = Orientation(p1, p2, q1);
        double o2 = Orientation(p1, p2, q2);
        double o3 = Orientation(q1, q2, p1);
        double o4 = Orientation(q1, q2, p2);

        // Standard intersection check
        if ((o1 * o2 < 0 && o3 * o4 < 0) || 
            (o1 == 0 && IsPointOnSegment(p1, p2, q2)) || 
            (o2 == 0 && IsPointOnSegment(p1, p2, q1)) ||
            (o3 == 0 && IsPointOnSegment(q1, q2, p1)) || 
            (o4 == 0 && IsPointOnSegment(q1, q2, p2))) 
        {
            // Compute exact intersection point
            double a1 = p2[1] - p1[1];
            double b1 = p1[0] - p2[0];
            double c1 = a1 * p1[0] + b1 * p1[1];

            double a2 = q2[1] - q1[1];
            double b2 = q1[0] - q2[0];
            double c2 = a2 * q1[0] + b2 * q1[1];

            double det = a1 * b2 - a2 * b1;
            
            if (std::abs(det) > epsilon) { // Avoid division by zero
                line_knot_span_edge_intersection_point[0] = (b2 * c1 - b1 * c2) / det;
                line_knot_span_edge_intersection_point[1] = (a1 * c2 - a2 * c1) / det;
                return true;
            }

            // Special case: return true if collinear and overlapping
            return true;
        }

        return false;
    }

    Vector ClassifyIntegrationPointsExtendedGradientMethodProcess::ComputeUnitNormal(
        const array_1d<double,3>& line_p1, const array_1d<double,3>& line_p2){
        Vector normal(3);
        normal[0] = - (line_p2[1] - line_p1[1]);
        normal[1] = (line_p2[0] - line_p1[0]);
        normal[2] = 0.0;
        return - normal/norm_2(normal);
    }

    std::pair<bool, int>  ClassifyIntegrationPointsExtendedGradientMethodProcess::IsRectangleInsidePolygon(const array_1d<double, 3>& p1, const array_1d<double, 3>& p2,
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

    // This method checks if a point is inside a polygon using the ray tracing algorithm
    bool ClassifyIntegrationPointsExtendedGradientMethodProcess::IsPointInsidePolygon(const array_1d<double, 3>& point, const std::vector<array_1d<double, 3>>& polygon_vertices) {
        int n = polygon_vertices.size();
        bool inside = false;

        // Loop over each edge of the polygon
        for (int i = 0; i < n; i++) {
            int j = (i == 0) ? (n - 1) : (i - 1);

            const array_1d<double, 3>& p1 = polygon_vertices[i];
            const array_1d<double, 3>& p2 = polygon_vertices[j];

            // Check if the point lies exactly on the edge
            if (IsPointOnSegment(point, p1, p2)) {
                return true; // The point is on the boundary
            }

            // Check if the ray from point crosses the edge from p1 to p2
            if (((p1[1] > point[1]) != (p2[1] > point[1])) && 
                (point[0] < (p2[0] - p1[0]) * (point[1] - p1[1]) / (p2[1] - p1[1]) + p1[0])) {
                inside = !inside;  // Toggle the inside flag
            }
        }
        return inside;  // Return whether the point is inside the polygon
    }

    // Helper function to check if a point lies on a line segment
    bool ClassifyIntegrationPointsExtendedGradientMethodProcess::IsPointOnSegment(const array_1d<double, 3>& point, 
    const array_1d<double, 3>& p1, const array_1d<double, 3>& p2){
         double cross_product = (point[1] - p1[1]) * (p2[0] - p1[0]) - (point[0] - p1[0]) * (p2[1] - p1[1]);

        // Check if the point is collinear with the edge (cross product ≈ 0)
        if (std::abs(cross_product) > 1e-9) {
            return false;
        }

        // Check if the point lies within the segment bounds
        double dot_product = (point[0] - p1[0]) * (p2[0] - p1[0]) + (point[1] - p1[1]) * (p2[1] - p1[1]);
        if (dot_product < 0) {
            return false;
        }

        double squared_length = (p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]);
        return dot_product <= squared_length;
    }

    // This method creates integration points over the inner and outer loops
    void ClassifyIntegrationPointsExtendedGradientMethodProcess::CreateBrepCurveOnSurfaceIntegrationPoints(std::shared_ptr<NurbsSurfaceGeometry<3, PointerVector<Node>>> p_nurbs_surface, 
    ModelPart& intersected_elements_sub_model_part, ModelPart* p_loop){
        // Create the integration points for the outer loop
        for (auto condition_it = p_loop->ConditionsBegin(); condition_it != p_loop->ConditionsEnd(); condition_it++){
            GeometryType condition_geometry = condition_it->GetGeometry();
            
            Geometry<GeometricalObject::NodeType>::PointsArrayType brep_curve_points;
            brep_curve_points.push_back(p_loop->pGetNode(condition_geometry.GetPoint(0).Id()));
            brep_curve_points.push_back(p_loop->pGetNode(condition_geometry.GetPoint(1).Id()));
                                
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
            IntegrationPointsArrayType integration_points;;
            IntegrationInfo integration_info = nurbs_curve_on_surface.GetDefaultIntegrationInfo();
            nurbs_curve_on_surface.CreateIntegrationPoints(integration_points, integration_info);
            
            // Create quadrature point geometries
            GeometriesArrayType quadrature_points_geometries(integration_points.size()); 
            IndexType number_of_shape_functions_derivatives = 2;
            nurbs_curve_on_surface.CreateQuadraturePointGeometries(quadrature_points_geometries, number_of_shape_functions_derivatives, integration_points, integration_info);
            
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
    }

} // End namespace Kratos
