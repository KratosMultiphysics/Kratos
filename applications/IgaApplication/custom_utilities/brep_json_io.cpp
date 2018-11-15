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
#include "iga_application.h"
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
            file << separator << "[" << element->Id() << ", " << element->GetValue(FACE_BREP_ID) << ",[" << local_parameters[0] << ", " << local_parameters[1] << "]]";
            separator = ",";
        }
        file << "]}";
        file.close();
    }

    std::vector<BrepModel> BrepJsonIO::ImportNurbsBrepGeometry(
        ModelPart& model_part,
        Parameters& rNurbsBrepGeometryJson)
    {
        std::cout << "\n> Start reading CAD geometry" << std::endl;

        double model_tolerance = rNurbsBrepGeometryJson["tolerances"]["model_tolerance"].GetDouble();

        std::vector<BrepModel> r_brep_model_vector;

        for (int brep_i = 0; brep_i < rNurbsBrepGeometryJson["breps"].size(); brep_i++)
        {
            Parameters brep_json = rNurbsBrepGeometryJson["breps"][brep_i];
            int brep_brep_id = brep_json["brep_id"].GetInt();

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

                //model_part.CreateSubModelPart("FACE_" + std::to_string(face_id));
                ModelPart& sub_model_part_face = model_part.CreateSubModelPart("FACE_" + std::to_string(face_id));

                //sub_model_part_face.CreateSubModelPart("FACE_" + std::to_string(face_id) + "_CPS");
                ModelPart& sub_model_part_face_cp = sub_model_part_face.CreateSubModelPart("FACE_" + std::to_string(face_id) + "_CPS");

                //std::cout << "> Reading face " << face_id << "..." << std::endl;

                bool is_trimmed = brep_json["faces"][i]["surface"]["is_trimmed"].GetBool();
                bool is_rational = brep_json["faces"][i]["surface"]["is_rational"].GetBool();

                // Variables needed later
                int length_u_vector = brep_json["faces"][i]["surface"]["knot_vectors"][0].size();
                int length_v_vector = brep_json["faces"][i]["surface"]["knot_vectors"][1].size();
                Vector knot_vector_u = ZeroVector(length_u_vector);
                Vector knot_vector_v = ZeroVector(length_v_vector);
                int p;
                int q;
                IntVector control_points_ids;

                // read and store knot_vector_u
                for (int u_idx = 0; u_idx < length_u_vector; u_idx++)
                {
                    knot_vector_u(u_idx) = brep_json["faces"][i]["surface"]["knot_vectors"][0][u_idx].GetDouble();
                }

                // read and store knot_vector_v
                for (int v_idx = 0; v_idx < length_v_vector; v_idx++)
                {
                    knot_vector_v(v_idx) = brep_json["faces"][i]["surface"]["knot_vectors"][1][v_idx].GetDouble();
                }

                // read and store polynamial degree p and q
                p = brep_json["faces"][i]["surface"]["degrees"][0].GetInt();
                q = brep_json["faces"][i]["surface"]["degrees"][1].GetInt();

                std::cout << "> Reading face " << face_id << " geometry data" << std::endl;

                // read and store control_points
                // Control points in each patch get a global as well as a mapping matrix id
                // brep Id: Unique Id for each control point (given by json-file)
                // mapping matrix id: specifies position in global mapping matrix (given by numbering of control points)
                for (int cp_idx = 0; cp_idx < brep_json["faces"][i]["surface"]["control_points"].size(); cp_idx++)
                {
                    unsigned int cp_id = brep_json["faces"][i]["surface"]["control_points"][cp_idx][0].GetInt();
                    double x = brep_json["faces"][i]["surface"]["control_points"][cp_idx][1][0].GetDouble();
                    double y = brep_json["faces"][i]["surface"]["control_points"][cp_idx][1][1].GetDouble();
                    double z = brep_json["faces"][i]["surface"]["control_points"][cp_idx][1][2].GetDouble();
                    double w = brep_json["faces"][i]["surface"]["control_points"][cp_idx][1][3].GetDouble();

                    control_points_ids.push_back(cp_id);

                    sub_model_part_face_cp.CreateNewNode(cp_id, x, y, z);
                    sub_model_part_face_cp.pGetNode(cp_id)->SetValue(NURBS_CONTROL_POINT_WEIGHT, w);
                }
                std::cout << "> Reading face " << face_id << " cps" << std::endl;

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // 2. step: boundary loops
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                std::vector<BrepBoundaryLoop> trimming_loops;
                //TrimmingCurveVector trimming_curves;

                // For better reading
                Parameters boundary_dict(brep_json["faces"][i]["boundary_loops"]);


                std::cout << "> Reading face " << face_id << " boundary loops" << std::endl;

                // Loop over all boundary loops
                for (int loop_idx = 0; loop_idx < boundary_dict.size(); loop_idx++)
                {
                    std::vector<BrepTrimmingCurve> loop_curves;

                    std::cout << "> Reading face " << face_id << " finishing" << std::endl;
                    // Loop over all curves
                    for (int edge_idx = 0; edge_idx < boundary_dict[loop_idx]["trimming_curves"].size(); edge_idx++)
                    {
                        int curve_index = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["trim_index"].GetInt();
                        //loop.push_back(curve_index);

                        std::cout << "> Reading face " << face_id << " and " << curve_index << " curve index." << std::endl;
                        bool curve_direction = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["curve_direction"].GetBool();
                        // Variables needed later
                        int length_u_vector = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["knot_vector"].size();
                        Vector boundary_knot_vector = ZeroVector(length_u_vector);
                        std::vector<BoundedVector<double, 4>> boundary_control_points;
                        Vector active_range = ZeroVector(2);
                        // read and store knot_vector_u
                        std::cout << "> Reading knot vector " << face_id << " " << std::endl;
                        for (int u_idx = 0; u_idx < length_u_vector; u_idx++)
                        {
                            boundary_knot_vector(u_idx) = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["knot_vector"][u_idx].GetDouble();
                        }

                        // read and store polynamial degree p and q
                        int boundary_p = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["degree"].GetInt();

                        std::cout << "> Reading control points " << face_id << " finishing" << std::endl;
                        // read and store control_points
                        for (int cp_idx = 0; cp_idx < boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"].size(); cp_idx++)
                        {
                            //unsigned int cp_id = extractInt(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][0]);
                            BoundedVector<double, 4> control_point;
                            control_point[0] = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1][0].GetDouble();
                            control_point[1] = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1][1].GetDouble();
                            control_point[2] = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1][2].GetDouble();
                            control_point[3] = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1][3].GetDouble();

                            boundary_control_points.push_back(control_point);
                        }

                        std::cout << "> Reading active range " << face_id << " finishing" << std::endl;
                        active_range(0) = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["active_range"][0].GetDouble();
                        active_range(1) = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["active_range"][1].GetDouble();

                        // Create and store edge
                        BrepTrimmingCurve new_boundary_curve(curve_index, boundary_knot_vector, boundary_p, boundary_control_points, curve_direction, true, active_range);

                        loop_curves.push_back(new_boundary_curve);
                    }
                    bool is_outer_loop = true;
                    if (boundary_dict[loop_idx]["loop_type"].GetString() == "inner")
                        is_outer_loop = false;
                    BrepBoundaryLoop loop(loop_curves, is_outer_loop);
                    trimming_loops.push_back(loop);
                }


                std::vector<BrepBoundaryLoop> embedded_loops;
                //TrimmingCurveVector trimming_curves;

                // For better reading
                Parameters embedded_dict(brep_json["faces"][i]["embedded_loops"]);

                // Loop over all embedded loops
                for (int loop_idx = 0; loop_idx < embedded_dict.size(); loop_idx++)
                {
                    std::vector<BrepTrimmingCurve> loop_curves;

                    // Loop over all curves
                    for (int edge_idx = 0; edge_idx < embedded_dict[loop_idx]["trimming_curves"].size(); edge_idx++)
                    {
                        int curve_index = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["trim_index"].GetInt();
                        //loop.push_back(curve_index);

                        bool curve_direction = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["curve_direction"].GetBool();
                        // Variables needed later
                        int length_u_vector = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["knot_vector"].size();
                        Vector boundary_knot_vector = ZeroVector(length_u_vector + 2);

                        std::vector<BoundedVector<double, 4>> boundary_control_points;
                        Vector active_range = ZeroVector(2);
                        // read and store knot_vector_u
                        for (int u_idx = 0; u_idx < length_u_vector; u_idx++)
                        {
                            boundary_knot_vector(u_idx) = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["knot_vector"][u_idx].GetDouble();
                        }

                        // read and store polynamial degree p and q
                        int boundary_p = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["degree"].GetInt();

                        // read and store control_points
                        for (int cp_idx = 0; cp_idx < embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"].size(); cp_idx++)
                        {
                            //unsigned int cp_id = extractInt(embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][0]);
                            BoundedVector<double, 4> control_point;
                            control_point[0] = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1][0].GetDouble();
                            control_point[1] = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1][1].GetDouble();
                            control_point[2] = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1][2].GetDouble();
                            control_point[3] = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1][3].GetDouble();

                            boundary_control_points.push_back(control_point);
                        }
                        active_range(0) = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["active_range"][0].GetDouble();
                        active_range(1) = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["active_range"][1].GetDouble();

                        // Create and store edge
                        BrepTrimmingCurve new_boundary_curve(curve_index, boundary_knot_vector, boundary_p, boundary_control_points, curve_direction, true, active_range);

                        loop_curves.push_back(new_boundary_curve);
                    }
                    bool is_outer_loop = true;
                    if (embedded_dict[loop_idx]["loop_type"].GetString() == "inner")
                        is_outer_loop = false;

                    std::cout << is_outer_loop << "length: "<< loop_curves.size() << std::endl;
                    BrepBoundaryLoop loop(loop_curves, is_outer_loop);
                    embedded_loops.push_back(loop);
                }

                std::vector<BrepFace::EmbeddedPoint> embedded_points;
                if (brep_json["faces"][i].Has("embedded_points"))
                {
                    // For better reading
                    Parameters embedded_points_dict(brep_json["faces"][i]["embedded_points"]);
                    for (int i = 0; i < embedded_points_dict.size(); i++)
                    {
                        int trim_index = embedded_points_dict[i]["trim_index"].GetInt();

                        Vector coords = ZeroVector(3);
                        coords(0) = embedded_points_dict[i]["point"][0].GetDouble();
                        coords(1) = embedded_points_dict[i]["point"][1].GetDouble();
                        coords(2) = embedded_points_dict[i]["point"][2].GetDouble();

                        BrepFace::EmbeddedPoint point(trim_index, coords);
                        embedded_points.push_back(point);
                    }
                }

                std::cout << "> Reading face " << face_id << " finishing" << std::endl;
                //ModelPart sub_model_part_face = sub_model_part_face.GetSubModelPart("FACE_" + std::to_string(face_id) + "_CPS");
                // create face
                BrepFace face(face_id, is_trimmed, is_rational, 
                    trimming_loops, embedded_loops, embedded_points, 
                    knot_vector_u, knot_vector_v, p, q,
                    control_points_ids, model_part);
                faces_vector.push_back(face);
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Edges
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            for (int i = 0; i < brep_json["edges"].size(); i++)
            {
                // For better reading
                Parameters edge_dict(brep_json["edges"][i]);

                std::vector<BrepEdge::EdgeTopology> brep_edge_topology_vector;

                int edge_id = edge_dict["brep_id"].GetInt();

                //model_part.CreateSubModelPart("FACE_" + std::to_string(face_id));
                ModelPart& sub_model_part_edge = model_part.CreateSubModelPart("EDGE_" + std::to_string(edge_id));

                std::cout << "> Reading edge " << edge_id << "..." << std::endl;

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // 3d curve
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                int degree = edge_dict["3d_curve"]["degree"].GetInt();

                // read and store active_range
                Vector active_range(2);
                active_range(0) = edge_dict["3d_curve"]["active_range"][0].GetDouble();
                active_range(1) = edge_dict["3d_curve"]["active_range"][1].GetDouble();

                // read and store knot_vector
                double length_knot_vector = edge_dict["3d_curve"]["knot_vector"].size();
                Vector knot_vector = ZeroVector(length_knot_vector);
                for (int u_idx = 0; u_idx < length_knot_vector; u_idx++)
                {
                    knot_vector(u_idx) = edge_dict["3d_curve"]["knot_vector"][u_idx].GetDouble();
                }


                IntVector control_points_ids_edge;
                // read and store control_points
                // Control points in each patch get a global as well as a mapping matrix id
                // brep Id: Unique Id for each control point (given by json-file)
                // mapping matrix id: specifies position in global mapping matrix (given by numbering of control points)
                for (int cp_idx = 0; cp_idx < edge_dict["3d_curve"]["control_points"].size(); cp_idx++)
                {
                    int cp_id = edge_dict["3d_curve"]["control_points"][cp_idx][0].GetInt();
                    double x = edge_dict["3d_curve"]["control_points"][cp_idx][1][0].GetDouble();
                    double y = edge_dict["3d_curve"]["control_points"][cp_idx][1][1].GetDouble();
                    double z = edge_dict["3d_curve"]["control_points"][cp_idx][1][2].GetDouble();
                    double w = edge_dict["3d_curve"]["control_points"][cp_idx][1][3].GetDouble();

                    control_points_ids_edge.push_back(cp_id);
                    sub_model_part_edge.CreateNewNode(cp_id, x, y, z);
                    sub_model_part_edge.GetNode(cp_id).SetValue(NURBS_CONTROL_POINT_WEIGHT, w);
                }

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // trimming range
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                std::vector<BrepEdge::TrimmingRange> brep_trimming_range_vector;
                if (brep_json["edges"].Has("trimmming_ranges"))
                {
                    for (int t = 0; t < brep_json["edges"]["trimmming_ranges"].size(); t++)
                    {
                        int trim_index = brep_json["edges"]["trimmming_ranges"][t]["trim_index"].GetInt();
                        Vector range(2);
                        range(0) = brep_json["edges"]["trimmming_ranges"][t]["range"][0].GetDouble();
                        range(1) = brep_json["edges"]["trimmming_ranges"][t]["range"][1].GetDouble();
                        BrepEdge::TrimmingRange trimming_range(trim_index, range);
                        brep_trimming_range_vector.push_back(trimming_range);
                    }
                }

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // edge topology
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                for (int j = 0; j < edge_dict["topology"].size(); j++)
                {
                    int face_id = edge_dict["topology"][j]["face_id"].GetInt();
                    int trim_index = edge_dict["topology"][j]["trim_index"].GetInt();
                    bool relative_direction = edge_dict["topology"][j]["relative_direction"].GetBool();

                    BrepEdge::EdgeTopology trim(face_id, trim_index, relative_direction);
                    brep_edge_topology_vector.push_back(trim);
                }

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // embedded points
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                std::vector<BrepEdge::EmbeddedPoint> brep_embedded_points_vector;
                if (edge_dict.Has("embedded_points"))
                {
                    for (int j = 0; j < edge_dict["embedded_points"].size(); j++)
                    {
                        int trim_index = edge_dict["embedded_points"][j]["trim_index"].GetInt();
                        double local_parameter = edge_dict["embedded_points"][j]["point"][0].GetDouble();

                        BrepEdge::EmbeddedPoint point(trim_index, local_parameter);
                        brep_embedded_points_vector.push_back(point);
                    }
                }

                std::cout << "> Finished reading edge " << edge_id << "..." << std::endl;
                BrepEdge edge(edge_id, brep_edge_topology_vector, brep_trimming_range_vector, brep_embedded_points_vector,
                    degree, knot_vector, active_range, control_points_ids_edge, model_part);

                std::cout << "> Finished creating edge " << edge_id << "..." << std::endl;
                edges_vector.push_back(edge);
            }
    //    //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    //    // Vertices
    //    //    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //    //    if (brep_json.Has("vertices"))
    //    //    {
    //    //        for (int i = 0; i < brep_json["vertices"].size(); i++)
    //    //        {
    //    //            // For better reading
    //    //            Parameters vertex_dict(brep_json["vertices"][i]);

    //    //            unsigned int vertex_id = vertex_dict["brep_id"].GetInt();

    //    //            std::cout << "> Start reading vertex " << vertex_id << "..." << std::endl;

    //    //            int cp_id = 0;// vertex_dict["coordinates"][0].GetInt();
    //    //            Vector coords = ZeroVector(4);
    //    //            //coords(0) = vertex_dict["coordinates"][1][0].GetDouble();
    //    //            //coords(1) = vertex_dict["coordinates"][1][1].GetDouble();
    //    //            //coords(2) = vertex_dict["coordinates"][1][2].GetDouble();
    //    //            //coords(3) = vertex_dict["coordinates"][1][3].GetDouble();

    //    //            std::vector<BrepVertex::Topology> brep_vertices_topology_vector;
    //    //            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    //            // topology
    //    //            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    //            for (int j = 0; j < vertex_dict["topology"].size(); j++)
    //    //            {
    //    //                unsigned int face_id = vertex_dict["topology"][j]["brep_id"].GetInt();
    //    //                unsigned int trim_index = vertex_dict["topology"][j]["trim_index"].GetInt();

    //    //                BrepVertex::Topology trim(face_id, trim_index);
    //    //                brep_vertices_topology_vector.push_back(trim);
    //    //            }
    //    //            std::cout << "> Finished reading vertex " << vertex_id << "..." << std::endl;

    //    //            BrepVertex vertex(vertex_id, brep_vertices_topology_vector, cp_id, coords);
    //    //            vertices_vector.push_back(vertex);
    //    //            std::cout << "> Finished creating vertex " << vertex_id << "..." << std::endl;
    //    //        }
    //    //    }
    //    //    std::cout << "Before creating BrepModel" << std::endl;

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // 3. Step: Create BrepModel
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            BrepModel brep(brep_brep_id, model_tolerance, faces_vector, edges_vector, vertices_vector);

            r_brep_model_vector.push_back(brep);// [brep_i] = &brep;
        }
        std::cout << "\n> Finished reading CAD geometry" << std::endl;
        return r_brep_model_vector;
    }

}  // namespace Kratos.
