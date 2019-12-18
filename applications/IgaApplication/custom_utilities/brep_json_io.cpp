//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "brep_json_io.h"
#include "iga_application_variables.h"


namespace Kratos
{
    void BrepJsonIO::WriteIntegrationDomainJson(ModelPart& rModelPart, const std::string& rOutputFileName)
    {
        std::ofstream file;
        file.open(rOutputFileName);
        std::string separator = "";
        file << "{\"geometry_integration_points\" :[";
        for (auto element = rModelPart.ElementsBegin(); element != rModelPart.ElementsEnd(); ++element)
        {
            const Vector& local_parameters = element->GetValue(LOCAL_COORDINATES);
            //const Vector& global_coordinates = element->Calculate(COORDINATES);
            file << separator << "[" << element->Id() << ", " << element->GetValue(BREP_ID) << ",[" << local_parameters[0] << ", " << local_parameters[1] << "]]";
            separator = ",";
        }
        file << "]}";
        file.close();
    }

    void BrepJsonIO::ExportNurbsGeometry(std::vector<BrepModel> m_brep_model_vector)
    {
        std::cout << "\n> Start writing CAD geometry" << std::endl;
        Parameters model_para; 
        
        // Model Tolerance
        Parameters model_tolerance_para; 
        const double model_tolerance =  m_brep_model_vector[0].GetModelTolerance();

        model_tolerance_para.AddEmptyValue("model_tolerance"); 
        model_tolerance_para["model_tolerance"].SetDouble(model_tolerance); 
        model_para.AddValue("tolerances", model_tolerance_para); 

        // Version Number
        model_para.AddEmptyValue("version_number");
        model_para["version_number"].SetInt(1); 

        // Loop over Breps
        for (int brep_i = 0; brep_i < m_brep_model_vector.size(); ++brep_i)
        {
            Parameters brep_para; 
            const int brep_id = m_brep_model_vector[brep_i].GetId();
            
            brep_para.AddEmptyValue("brep_id");
            brep_para["brep_id"].SetInt(brep_id); 

            // Loop over Faces
            for (int face_i = 0; face_i < m_brep_model_vector[brep_i].GetFaceVector().size(); ++face_i)
            {
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // 1. Step: faces
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                
                Parameters faces_para; 
                const int face_id = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetId(); 
                const bool swapped_surface_normal = false; 

                faces_para.AddEmptyValue("brep_id");
                faces_para.AddEmptyValue("swapped_surface_normal"); 
                faces_para["brep_id"].SetInt(face_id);
                faces_para["swapped_surface_normal"].SetBool(swapped_surface_normal);  
                
                Parameters surface_para; 
                const bool is_trimmed = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetIsTrimmed(); 
                const bool is_rational = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetIsRational();     

                Vector degrees = ZeroVector(2);
                degrees[0] = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->DegreeU();
                degrees[1] = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->DegreeV();
                
                const std::vector<double> knots_u = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->KnotsU(); 
                const std::vector<double> knots_v = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->KnotsV();

                const int number_knots_u = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->NbKnotsU();
                const int number_knots_v = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->NbKnotsV();

                Vector knot_vector_u = ZeroVector(number_knots_u + 2); 
                Vector knot_vector_v = ZeroVector(number_knots_v + 2); 

                knot_vector_u[0] = knots_u[0];
                knot_vector_v[0] = knots_v[0];
                for(int i = 0; i < number_knots_u; ++i)    knot_vector_u[i+1] = knots_u[i]; 
                for(int i = 0; i < number_knots_v; ++i)    knot_vector_v[i+1] = knots_v[i]; 
                knot_vector_u[number_knots_u + 1] = knots_u[number_knots_u - 1];
                knot_vector_v[number_knots_v + 1] = knots_v[number_knots_v - 1]; 

                const int number_cps_u = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->NbKnotsU();
                const int number_cps_v = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->NbKnotsV();

                surface_para.AddEmptyValue("is_trimmed"); 
                surface_para.AddEmptyValue("is_rational"); 
                surface_para.AddEmptyArray("degrees");
                surface_para.AddEmptyArray("knot_vectors");
                surface_para.AddEmptyArray("control_points");
                surface_para["is_trimmed"].SetBool(is_trimmed); 
                surface_para["is_rational"].SetBool(is_rational); 
                surface_para["degrees"].SetVector(degrees); 
                surface_para["knot_vectors"].Append(knot_vector_u); 
                surface_para["knot_vectors"].Append(knot_vector_v); 
                
                Vector control_points = ZeroVector(4); 
                for (int i = 0; i < number_cps_u; ++i)
                {
                    for (int j = 0; j < number_cps_v; ++j)
                    {
                        const int control_points_ids = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->GetNode(i, j)->GetId(); 
                        control_points[0] = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->GetNode(j,i)->X(); 
                        control_points[1] = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->GetNode(j,i)->Y(); 
                        control_points[2] = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->GetNode(j,i)->Z(); 
                        control_points[3] = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetSurface()->GetNode(j,i)->GetValue(NURBS_CONTROL_POINT_WEIGHT); 
                        
                        Parameters control_points_para; 
                        control_points_para.AddEmptyArray("cp"); 
                        control_points_para["cp"].Append(control_points_ids); 
                        control_points_para["cp"].Append(control_points); 
                        surface_para["control_points"].Append(control_points_para["cp"]); 
                    }
                }

                faces_para.AddValue("surface", surface_para); 

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // 2. step: boundary loops
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                auto boundary_loop = m_brep_model_vector[brep_i].GetFaceVector()[face_i].GetBoundaryLoop(); 
                
                // loop over boundaries
                for (int b_loop_i = 0; b_loop_i < boundary_loop.size(); ++b_loop_i)
                {
                    Parameters boundary_loop_para; 
                    boundary_loop_para.AddEmptyValue("loop_type"); 

                    if (boundary_loop[b_loop_i].IsOuterLoop() == true)
                    {
                        boundary_loop_para["loop_type"].SetString("outer"); 
                    }   
                    else
                    {
                        boundary_loop_para["loop_type"].SetString("inner"); 
                    }   

                    auto trimming_curves = boundary_loop[b_loop_i].GetTrimmingCurves();
                    
                    for (int tcurves_i = 0; tcurves_i < trimming_curves.size(); ++tcurves_i)
                    {
                        Parameters trimming_curves_para;
                        trimming_curves_para.AddEmptyValue("trim_index");
                        trimming_curves_para.AddEmptyValue("curve_direction");
                        trimming_curves_para["trim_index"].SetInt(trimming_curves[tcurves_i].GetTrimIndex());
                        trimming_curves_para["curve_direction"].SetBool(trimming_curves[tcurves_i].GetCurveDirection());

                        Parameters parameter_curve_para;

                        const int parameter_curve_degree = trimming_curves[tcurves_i].GetCurve2D()->Degree();
                        const bool parameter_curve_isrational = trimming_curves[tcurves_i].GetCurve2D()->CurveGeometry()->IsRational(); 
                        const std::vector<double> parameter_curve_knots = trimming_curves[tcurves_i].GetCurve2D()->CurveGeometry()->Knots(); 
                        const int parameter_curve_nb_knots = trimming_curves[tcurves_i].GetCurve2D()->CurveGeometry()->NbKnots(); 

                        Vector parameter_curve_knot_vector = ZeroVector(parameter_curve_nb_knots + 2); 

                        parameter_curve_knot_vector[0] = parameter_curve_knots[0];
                        for(int i = 0; i < parameter_curve_nb_knots; ++i)    parameter_curve_knot_vector[i+1] = parameter_curve_knots[i];
                        parameter_curve_knot_vector[parameter_curve_nb_knots + 1] = parameter_curve_knots[parameter_curve_nb_knots - 1];

                        Vector boundary_active_range = ZeroVector(2); 
                        boundary_active_range[0] = parameter_curve_knots.front(); 
                        boundary_active_range[1] = parameter_curve_knots.back();  

                        parameter_curve_para.AddEmptyValue("is_rational"); 
                        parameter_curve_para.AddEmptyValue("degree"); 
                        parameter_curve_para.AddEmptyArray("knot_vector"); 
                        parameter_curve_para.AddEmptyArray("active_range");
                        parameter_curve_para.AddEmptyArray("control_points");
                        
                        parameter_curve_para["is_rational"].SetBool(parameter_curve_isrational);
                        parameter_curve_para["degree"].SetInt(parameter_curve_degree); 
                        parameter_curve_para["knot_vector"].SetVector(parameter_curve_knot_vector);
                        parameter_curve_para["active_range"].SetVector(boundary_active_range);
                        

                        const int number_cps = trimming_curves[tcurves_i].GetCurve2D()->CurveGeometry()->NbPoles(); 

                        Vector control_points = ZeroVector(4); 

                        for(int i = 0; i < number_cps; ++i)
                        {
                            control_points[0] = trimming_curves[tcurves_i].GetCurve2D()->CurveGeometry()->Poles()[i][0];
                            control_points[1] = trimming_curves[tcurves_i].GetCurve2D()->CurveGeometry()->Poles()[i][1];
                            control_points[2] = 0; 
                            control_points[3] = trimming_curves[tcurves_i].GetCurve2D()->CurveGeometry()->Weight(i); 
                            
                            Parameters control_points_para; 
                            control_points_para.AddEmptyArray("cp"); 
                            control_points_para["cp"].Append(0); 
                            control_points_para["cp"].Append(control_points); 
                            parameter_curve_para["control_points"].Append(control_points_para["cp"]); 
                        }
                        
                        trimming_curves_para.AddValue("parameter_curve", parameter_curve_para); 
                        boundary_loop_para.AddEmptyArray("trimming_curves"); 
                        boundary_loop_para["trimming_curves"].Append(trimming_curves_para);                         
                    }

                    faces_para.AddEmptyArray("boundary_loops"); 
                    faces_para["boundary_loops"].Append(boundary_loop_para);

                    Vector embedded_loops = ZeroVector(0);
                    Vector embedded_points = ZeroVector(0); 

                    faces_para.AddEmptyValue("embedded_loops"); 
                    faces_para.AddEmptyValue("embedded_points"); 
                    
                    faces_para["embedded_loops"].SetVector(embedded_loops);
                    faces_para["embedded_points"].SetVector(embedded_points);
                } 

                brep_para.AddEmptyArray("faces"); 
                brep_para["faces"].Append(faces_para); 
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Edges
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            //loop over edges
            for (int edge_i = 0; edge_i < m_brep_model_vector[brep_i].GetEdgeVector().size(); ++edge_i)
            {
                Parameters edges_para; 
                const int edge_id = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetId(); 

                edges_para.AddEmptyValue("brep_id"); 
                edges_para["brep_id"].SetInt(edge_id); 
                
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // 3d curve
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                
                Parameters curve_3d_para;

                const int degree = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->Degree(); 
                const std::vector<double> knots = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->Knots(); 
                const int number_knots = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->NbKnots();
                
                Vector knot_vector = ZeroVector(number_knots + 2); 

                knot_vector[0] = knots[0];
                for(int i = 0; i < number_knots; ++i)    knot_vector[i+1] = knots[i];
                knot_vector[number_knots + 1] = knots[number_knots - 1];

                Vector curve_3d_active_range = ZeroVector(2);
                curve_3d_active_range[0] = knots.front(); 
                curve_3d_active_range[1] = knots.back(); 

                const int number_cps = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->NbPoles(); 
                
                curve_3d_para.AddEmptyValue("degree");
                curve_3d_para.AddEmptyArray("knot_vector");
                curve_3d_para.AddEmptyArray("active_range");
                curve_3d_para.AddEmptyArray("control_points");
                curve_3d_para["degree"].SetInt(degree);
                curve_3d_para["knot_vector"].SetVector(knot_vector);
                curve_3d_para["active_range"].SetVector(curve_3d_active_range);
            
                Vector control_points = ZeroVector(4); 
            
                for (int i = 0; i < number_cps; ++i)
                {
                    const int control_points_ids = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->GetNode(i)->GetId(); 
                    control_points[0] = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->GetNode(i)->X(); 
                    control_points[1] = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->GetNode(i)->Y(); 
                    control_points[2] = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->GetNode(i)->Z(); 
                    control_points[3] = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetCurve3d()->GetNode(i)->GetValue(NURBS_CONTROL_POINT_WEIGHT); 
                
                    Parameters control_points_para;     
                    control_points_para.AddEmptyArray("cp"); 
                    control_points_para["cp"].Append(control_points_ids); 
                    control_points_para["cp"].Append(control_points); 
                    curve_3d_para["control_points"].Append(control_points_para["cp"]); 
                }
                edges_para.AddValue("3d_curve", curve_3d_para);

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // edge topology
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                Parameters topology_para;
                
                for (int i = 0; i < m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetBrepEdgeTopologyVector().size(); ++i)
                {
                    topology_para.AddEmptyValue("face_id"); 
                    topology_para.AddEmptyValue("trim_index");
                    topology_para.AddEmptyValue("relative_direction");
                    const int topology_face_id = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetBrepEdgeTopologyVector()[i].brep_id;
                    const int topology_trim_index = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetBrepEdgeTopologyVector()[i].trim_index;
                    const bool topology_relative_direction = m_brep_model_vector[brep_i].GetEdgeVector()[edge_i].GetBrepEdgeTopologyVector()[i].relative_direction;
                    
                    topology_para["face_id"].SetInt(topology_face_id); 
                    topology_para["trim_index"].SetInt(topology_trim_index); 
                    topology_para["relative_direction"].SetBool(topology_relative_direction); 
                }

                edges_para.AddEmptyArray("topology"); 
                edges_para["topology"].Append(topology_para);

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // embedded points
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                edges_para.AddEmptyArray("embedded_points"); 
                brep_para.AddEmptyArray("edges"); 
                brep_para["edges"].Append(edges_para); 
            }
            
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Vertices
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            brep_para.AddEmptyArray("vertices"); 
            model_para.AddEmptyArray("breps");
            model_para["breps"].Append(brep_para); 
        }

        std::cout << model_para << std::endl;

    }

    std::vector<BrepModel> BrepJsonIO::ImportNurbsBrepGeometry(
        ModelPart& rModelPart,
        Parameters rNurbsBrepGeometryJson)
    {
        KRATOS_INFO("IGA") << "Start import CAD geometries" << std::endl;

        const double model_tolerance = rNurbsBrepGeometryJson["tolerances"]["model_tolerance"].GetDouble();
        std::vector<BrepModel> brep_model_vector;

        for (int brep_i = 0; brep_i < rNurbsBrepGeometryJson["breps"].size(); brep_i++)
        {
            const Parameters& brep_json = rNurbsBrepGeometryJson["breps"][brep_i];
            const int brep_brep_id = brep_json["brep_id"].GetInt();

            std::vector<BrepFace>   faces_vector;
            std::vector<BrepEdge>   edges_vector;
            std::vector<BrepVertex> vertices_vector;

            // loop over faces
            for (int i = 0; i < brep_json["faces"].size(); i++)
            {

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // 1. Step: faces
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                int face_id = brep_json["faces"][i]["brep_id"].GetInt();
                KRATOS_INFO_IF("IGA", mEchoLevel >= 1) << "Reading face " << face_id << "..." << std::endl;


                //model_part.CreateSubModelPart("FACE_" + std::to_string(face_id));
                ModelPart& sub_model_part_face = rModelPart.CreateSubModelPart("FACE_" + std::to_string(face_id));

                //sub_model_part_face.CreateSubModelPart("FACE_" + std::to_string(face_id) + "_CPS");
                ModelPart& sub_model_part_face_cp = sub_model_part_face.CreateSubModelPart("FACE_" + std::to_string(face_id) + "_CPS");

                //std::cout << "> Reading face " << face_id << "..." << std::endl;

                bool is_trimmed = brep_json["faces"][i]["surface"]["is_trimmed"].GetBool();
                bool is_rational = brep_json["faces"][i]["surface"]["is_rational"].GetBool();

                Vector knot_vector_u = brep_json["faces"][i]["surface"]["knot_vectors"][0].GetVector();
                Vector knot_vector_v = brep_json["faces"][i]["surface"]["knot_vectors"][1].GetVector();

                int p = brep_json["faces"][i]["surface"]["degrees"][0].GetInt();
                int q = brep_json["faces"][i]["surface"]["degrees"][1].GetInt();

                IntVector control_points_ids;

                KRATOS_INFO_IF("IGA", mEchoLevel >= 2) << "Reading face " << face_id << " cps" << std::endl;
                // read and store control_points
                // Control points in each patch get a global as well as a mapping matrix id
                // brep Id: Unique Id for each control point (given by json-file)
                // mapping matrix id: specifies position in global mapping matrix (given by numbering of control points)
                for (int cp_idx = 0; cp_idx < brep_json["faces"][i]["surface"]["control_points"].size(); cp_idx++)
                {
                    unsigned int cp_id = brep_json["faces"][i]["surface"]["control_points"][cp_idx][0].GetInt();
                    Vector cp = brep_json["faces"][i]["surface"]["control_points"][cp_idx][1].GetVector();

                    control_points_ids.push_back(cp_id);

                    sub_model_part_face_cp.CreateNewNode(cp_id, cp[0], cp[1], cp[2]);
                    sub_model_part_face_cp.pGetNode(cp_id)->SetValue(NURBS_CONTROL_POINT_WEIGHT, cp[3]);
                }
                KRATOS_INFO_IF("IGA", mEchoLevel >= 3) << "Reading face " << face_id << " cps finished" << std::endl;

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // 2. step: boundary loops
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                std::vector<BrepBoundaryLoop> trimming_loops;
                //TrimmingCurveVector trimming_curves;

                // For better reading
                const Parameters& boundary_dict(brep_json["faces"][i]["boundary_loops"]);


                KRATOS_INFO_IF("IGA", mEchoLevel >= 2) << "Reading face " << face_id << " boundary loops" << std::endl;
                // Loop over all boundary loops
                for (int loop_idx = 0; loop_idx < boundary_dict.size(); loop_idx++)
                {
                    std::vector<BrepTrimmingCurve> loop_curves;

                    // Loop over all curves
                    for (int trim_idx = 0; trim_idx < boundary_dict[loop_idx]["trimming_curves"].size(); trim_idx++)
                    {
                        ImportTrimmingCurve(boundary_dict[loop_idx]["trimming_curves"][trim_idx], loop_curves);
                    }
                    bool is_outer_loop = true;
                    if (boundary_dict[loop_idx]["loop_type"].GetString() == "inner")
                        is_outer_loop = false;
                    BrepBoundaryLoop loop(loop_curves, is_outer_loop);
                    trimming_loops.push_back(loop);
                }
                KRATOS_INFO_IF("IGA", mEchoLevel >= 3) << "Reading face " << face_id << " boundary loops finished" << std::endl;


                std::vector<BrepBoundaryLoop> embedded_loops;
                //TrimmingCurveVector trimming_curves;

                Parameters embedded_dict(brep_json["faces"][i]["embedded_loops"]);

                // Loop over all embedded loops
                for (std::size_t loop_idx = 0; loop_idx < embedded_dict.size(); loop_idx++)
                {
                    std::vector<BrepTrimmingCurve> loop_curves;

                    // Loop over all curves
                    for (std::size_t curve_idx = 0; curve_idx < embedded_dict[loop_idx]["trimming_curves"].size(); curve_idx++)
                    {
                        ImportTrimmingCurve(embedded_dict[loop_idx]["trimming_curves"][curve_idx], loop_curves);
                    }
                    bool is_outer_loop = true;
                    if (embedded_dict[loop_idx]["loop_type"].GetString() == "inner")
                        is_outer_loop = false;

                    std::cout << is_outer_loop << "length: " << loop_curves.size() << std::endl;
                    BrepBoundaryLoop loop(loop_curves, is_outer_loop);
                    embedded_loops.push_back(loop);
                }

                std::vector<BrepTrimmingCurve> trimming_curves;
                if (brep_json["faces"][i].Has("embedded_edges"))
                {
                    const Parameters& embedded_edges_dict(brep_json["faces"][i]["embedded_edges"]);
                    for (int i = 0; i < embedded_edges_dict.size(); ++i)
                    {
                        ImportTrimmingCurve(embedded_edges_dict[i], trimming_curves);
                    }
                }

                

                std::vector<BrepFace::EmbeddedPoint> embedded_points;
                if (brep_json["faces"][i].Has("embedded_points"))
                {
                    // For better reading
                    Parameters embedded_points_dict(brep_json["faces"][i]["embedded_points"]);
                    for (std::size_t i = 0; i < embedded_points_dict.size(); i++)
                    {
                        int trim_index = embedded_points_dict[i]["trim_index"].GetInt();

                        Vector coords = embedded_points_dict[i]["point"].GetVector();

                        BrepFace::EmbeddedPoint point(trim_index, coords);
                        embedded_points.push_back(point);
                    }
                }

                // create face
                BrepFace face(
                    face_id,
                    is_trimmed,
                    is_rational,
                    trimming_loops,
                    embedded_loops,
                    trimming_curves,
                    embedded_points,
                    knot_vector_u,
                    knot_vector_v,
                    p,
                    q,
                    control_points_ids,
                    rModelPart);
                faces_vector.push_back(face);
            }
            
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Edges
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (brep_json.Has("edges"))
            {
                ImportBrepEdges(
                    brep_json["edges"],
                    edges_vector,
                    rModelPart);
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Vertices
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            if (brep_json.Has("vertices"))
            {
                ImportBrepVertices(
                    brep_json["vertices"],
                    vertices_vector,
                    rModelPart);
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // 3. Step: Create BrepModel
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            BrepModel brep(
                brep_brep_id,
                model_tolerance,
                faces_vector,
                edges_vector,
                vertices_vector);
            brep_model_vector.push_back(brep);// [brep_i] = &brep;
        }
        KRATOS_INFO("IGA") << "Finished import CAD geometries" << std::endl;
        return brep_model_vector;
    }

    void BrepJsonIO::ImportBrepEdges(
        const Parameters& rEdges,
        std::vector<BrepEdge>& rEdgesVector,
        ModelPart& rModelPart)
    {
        for (std::size_t i = 0; i < rEdges.size(); i++)
        {
            Parameters edge_dict = rEdges[i];

            std::vector<BrepEdge::EdgeTopology> brep_edge_topology_vector;
            int edge_id = edge_dict["brep_id"].GetInt();
            // KRATOS_WATCH(edge_id)
            ModelPart& sub_model_part_edge = rModelPart.CreateSubModelPart("EDGE_" + std::to_string(edge_id));

            KRATOS_INFO_IF("IGA", mEchoLevel >= 2) << "Reading edge " << edge_id << "..." << std::endl;

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // 3d curve
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            int degree = edge_dict["3d_curve"]["degree"].GetInt();

            Vector active_range = edge_dict["3d_curve"]["active_range"].GetVector();

            // read and store knot_vector
            double length_knot_vector = edge_dict["3d_curve"]["knot_vector"].size();
            Vector knot_vector = edge_dict["3d_curve"]["knot_vector"].GetVector();

            IntVector control_points_ids_edge;
            // read and store control_points
            // Control points in each patch get a global as well as a mapping matrix id
            // brep Id: Unique Id for each control point (given by json-file)
            // mapping matrix id: specifies position in global mapping matrix (given by numbering of control points)
            for (std::size_t cp_idx = 0; cp_idx < edge_dict["3d_curve"]["control_points"].size(); cp_idx++)
            {
                int cp_id = edge_dict["3d_curve"]["control_points"][cp_idx][0].GetInt();
                Vector coordinates = edge_dict["3d_curve"]["control_points"][cp_idx][1].GetVector();

                control_points_ids_edge.push_back(cp_id);
                sub_model_part_edge.CreateNewNode(cp_id, coordinates[0], coordinates[1], coordinates[2]);
                sub_model_part_edge.GetNode(cp_id).SetValue(NURBS_CONTROL_POINT_WEIGHT, coordinates[3]);
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // trimming range
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            std::vector<BrepEdge::TrimmingRange> brep_trimming_range_vector;
            // if (brep_json["edges"].Has("trimmming_ranges"))
            // {
            //     for (std::size_t t = 0; t < brep_json["edges"]["trimmming_ranges"].size(); t++)
            //     {
            //         int trim_index = brep_json["edges"]["trimmming_ranges"][t]["trim_index"].GetInt();
            //         Vector range(2);
            //         range(0) = brep_json["edges"]["trimmming_ranges"][t]["range"][0].GetDouble();
            //         range(1) = brep_json["edges"]["trimmming_ranges"][t]["range"][1].GetDouble();
            //         BrepEdge::TrimmingRange trimming_range(trim_index, range);
            //         brep_trimming_range_vector.push_back(trimming_range);
            //     }
            // }
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // edge topology
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            if (edge_dict.Has("topology"))
            {
                for (std::size_t j = 0; j < edge_dict["topology"].size(); j++)
                {
                    int face_id = edge_dict["topology"][j]["brep_id"].GetInt();
                    int trim_index = edge_dict["topology"][j]["trim_index"].GetInt();
                    bool relative_direction = edge_dict["topology"][j]["relative_direction"].GetBool();

                    BrepEdge::EdgeTopology trim(face_id, trim_index, relative_direction);
                    brep_edge_topology_vector.push_back(trim);
                }
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // embedded points
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            std::vector<BrepEdge::EmbeddedPoint> brep_embedded_points_vector;
            if (edge_dict.Has("embedded_points"))
            {
                for (std::size_t j = 0; j < edge_dict["embedded_points"].size(); j++)
                {
                    int trim_index = edge_dict["embedded_points"][j]["trim_index"].GetInt();
                    double local_parameter = edge_dict["embedded_points"][j]["point"][0].GetDouble();

                    BrepEdge::EmbeddedPoint point(trim_index, local_parameter);
                    brep_embedded_points_vector.push_back(point);
                }
            }

            BrepEdge edge(
                edge_id,
                brep_edge_topology_vector,
                brep_trimming_range_vector,
                brep_embedded_points_vector,
                degree,
                knot_vector,
                active_range,
                control_points_ids_edge,
                rModelPart);

            rEdgesVector.push_back(edge);
        }
    }

    void BrepJsonIO::ImportBrepVertices(
        const Parameters& rVertices,
        std::vector<BrepVertex>& rVerticesVector,
        ModelPart& rModelPart)
    {
        for (int i = 0; i < rVertices.size(); i++)
        {
            Parameters vertex_dict(rVertices[i]);

            unsigned int vertex_id = vertex_dict["brep_id"].GetInt();

            KRATOS_INFO_IF("IGA", mEchoLevel >= 2) << "Reading vertex " << vertex_id << "..." << std::endl;

            int cp_id = 0;
            Vector coords = ZeroVector(3);

            if (vertex_dict.Has("coordinates"))
            {
                int cp_id = vertex_dict["coordinates"][0].GetInt();
                Vector coords = vertex_dict["coordinates"][1].GetVector();
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // topology
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            std::vector<BrepVertex::VertexTopology> brep_vertices_topology_vector;
            if (vertex_dict.Has("topology"))
            {
                for (int j = 0; j < vertex_dict["topology"].size(); j++)
                {
                    unsigned int brep_id = vertex_dict["topology"][j]["brep_id"].GetInt();
                    unsigned int trim_index = vertex_dict["topology"][j]["trim_index"].GetInt();

                    KRATOS_INFO_IF("IGA", mEchoLevel >= 4) << "Reading vertex topology of " << vertex_id 
                        << ". brep_id: " << brep_id << ", trim_index: " << trim_index << std::endl;

                    BrepVertex::VertexTopology trim(brep_id, trim_index);
                    brep_vertices_topology_vector.push_back(trim);
                }
            }

            BrepVertex vertex(vertex_id, brep_vertices_topology_vector, cp_id, coords);
            rVerticesVector.push_back(vertex);
        }
    }


    void BrepJsonIO::ImportTrimmingCurve(
        const Parameters& rTrimmingCurve,
        std::vector<BrepTrimmingCurve>& rTrimmingCurves)
    {
        int trim_index = rTrimmingCurve["trim_index"].GetInt();

        bool curve_direction = rTrimmingCurve["curve_direction"].GetBool();
        bool is_rational = rTrimmingCurve["parameter_curve"]["is_rational"].GetBool();

        Vector boundary_knot_vector = rTrimmingCurve["parameter_curve"]["knot_vector"].GetVector();

        // read and store polynamial degree p
        int degree = rTrimmingCurve["parameter_curve"]["degree"].GetInt();
        std::vector<BoundedVector<double, 4>> boundary_control_points;
        // read and store control_points
        for (std::size_t cp_idx = 0; cp_idx < rTrimmingCurve["parameter_curve"]["control_points"].size(); cp_idx++)
        {
            BoundedVector<double, 4> control_point = rTrimmingCurve["parameter_curve"]["control_points"][cp_idx][1].GetVector();
            boundary_control_points.push_back(control_point);
        }

        Vector active_range = rTrimmingCurve["parameter_curve"]["active_range"].GetVector();
        KRATOS_ERROR_IF(active_range.size() > 2) << "active_range of trim curve " << trim_index << " does not has the size 2" << std::endl;

        // Create and store edge
        BrepTrimmingCurve trimming_curve(
            trim_index,
            boundary_knot_vector,
            degree,
            boundary_control_points,
            curve_direction,
            is_rational,
            active_range);

        rTrimmingCurves.push_back(trimming_curve);
    }
}  // namespace Kratos.
