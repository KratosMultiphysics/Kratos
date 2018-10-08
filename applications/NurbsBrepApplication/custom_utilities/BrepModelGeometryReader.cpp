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
//#include "utilities/math_utils.h"


// External includes


// Project includes
#include "BrepModelGeometryReader.h"
#include "nurbs_brep_application.h"
#include "nurbs_brep_application_variables.h"


namespace Kratos
{
//TODO: you should isolate this from python. In kratos the "Parameters" class gives you access to the json object.
  //std::vector<BrepModel>
  void BrepModelGeometryReader::WriteGaussPoints(ModelPart& model_part)
  {
    ModelPart& faces = model_part.GetSubModelPart("FACES");
    std::ofstream file;
    file.open("faces.txt");
    for (auto node = faces.NodesBegin(); node != faces.NodesEnd(); ++node)
    {
      file << node->X() << " " << node->Y() << " " << node->Z() << " " << "\n";
    }
    file.close();

	if (faces.HasSubModelPart("FACE_2_EMBEDDED"))
	{
		ModelPart& faces_2_embedded = faces.GetSubModelPart("FACE_2_EMBEDDED");
		std::ofstream file_embedded;
		file_embedded.open("facesembedded.txt");
		for (auto node = faces_2_embedded.NodesBegin(); node != faces_2_embedded.NodesEnd(); ++node)
		{
			file_embedded << node->X() << " " << node->Y() << " " << node->Z() << " " << "\n";
		}
		file_embedded.close();
	}

	if (faces.HasSubModelPart("FACE_2_REVERSE"))
	{
		ModelPart& faces_2_reversed = faces.GetSubModelPart("FACE_2_REVERSE");
		std::ofstream file_reversed;
		file_reversed.open("facesreversed.txt");
		for (auto node = faces_2_reversed.NodesBegin(); node != faces_2_reversed.NodesEnd(); ++node)
		{
			file_reversed << node->X() << " " << node->Y() << " " << node->Z() << " " << "\n";
		}
		file_reversed.close();
	}

    ModelPart& edges = model_part.GetSubModelPart("EDGES");
    std::ofstream file_edges;
    file_edges.open("edges.txt");
    for (auto node = edges.NodesBegin(); node != edges.NodesEnd(); ++node)
    {
      file_edges << node->X() << " " << node->Y() << " " << node->Z() << " " << "\n";
    }
    file_edges.close();

    ModelPart& c_edges = model_part.GetSubModelPart("COUPLING_EDGES");
    std::ofstream file_c_edges;
    file_c_edges.open("c_edges.txt");
    for (auto node = c_edges.NodesBegin(); node != c_edges.NodesEnd(); ++node)
    {
      file_c_edges << node->X() << " " << node->Y() << " " << node->Z() << " ";
      
      Vector location = node->GetValue(LOCATION_SLAVE);
      file_c_edges << location[0] << " " << location[1] << " " << location[2] << "\n";
    }
    file_c_edges.close();
  }

  void BrepModelGeometryReader::WriteGaussPointsIteration(ModelPart& model_part, int step)
  {
      //ModelPart& faces = model_part.GetSubModelPart("FACES");
      std::ofstream file;
      file.open("faces" + std::to_string(step) + ".txt");
      for (auto node = model_part.ElementsBegin(); node != model_part.ElementsEnd(); ++node)
      {
          Vector variable = ZeroVector(3);
          node->Calculate(COORDINATES, variable, model_part.GetProcessInfo());
          file << node->Id() << " " << variable[0] << " " << variable[1] << " " << variable[2] << " " << "\n";
      }
      file.close();

      //if (faces.HasSubModelPart("FACE_2_EMBEDDED"))
      //{
      //    ModelPart& faces_2_embedded = faces.GetSubModelPart("FACE_2_EMBEDDED");
      //    std::ofstream file_embedded;
      //    file_embedded.open("facesembedded.txt");
      //    for (auto node = faces_2_embedded.NodesBegin(); node != faces_2_embedded.NodesEnd(); ++node)
      //    {
      //        file_embedded << node->X() << " " << node->Y() << " " << node->Z() << " " << "\n";
      //    }
      //    file_embedded.close();
      //}

      //if (faces.HasSubModelPart("FACE_2_REVERSE"))
      //{
      //    ModelPart& faces_2_reversed = faces.GetSubModelPart("FACE_2_REVERSE");
      //    std::ofstream file_reversed;
      //    file_reversed.open("facesreversed.txt");
      //    for (auto node = faces_2_reversed.NodesBegin(); node != faces_2_reversed.NodesEnd(); ++node)
      //    {
      //        file_reversed << node->X() << " " << node->Y() << " " << node->Z() << " " << "\n";
      //    }
      //    file_reversed.close();
      //}

      //const ModelPart& edges = model_part.GetSubModelPart("EDGES");
      std::ofstream file_edges;
      file_edges.open("edges" + std::to_string(step) + ".txt");
      for (auto node = model_part.ConditionsBegin(); node != model_part.ConditionsEnd(); ++node)
      {
          Vector variable = ZeroVector(3);
          node->Calculate(COORDINATES, variable, model_part.GetProcessInfo());
          file_edges << variable[0] << " " << variable[1] << " " << variable[2] << " " << "\n";
      }
      file_edges.close();

      //const ModelPart& c_edges = model_part.GetSubModelPart("COUPLING_EDGES");
      //std::ofstream file_c_edges;
      //file_c_edges.open("c_edges.txt");
      //for (auto node = c_edges.NodesBegin(); node != c_edges.NodesEnd(); ++node)
      //{
      //    file_c_edges << node->X() << " " << node->Y() << " " << node->Z() << " ";

      //    Vector location = node->GetValue(LOCATION_SLAVE);
      //    file_c_edges << location[0] << " " << location[1] << " " << location[2] << "\n";
      //}
      //file_c_edges.close();
  }

  void BrepModelGeometryReader::WriteGaussPointsJson(ModelPart& rModelPart, const std::string& rOutputFileName)
  {
	  ModelPart& faces = rModelPart.GetSubModelPart("FACES");
	  std::ofstream file;
	  file.open(rOutputFileName);
	  std::string separator = "";
	  file << "{\"gauss_points\" :[";
	  for (auto node = faces.NodesBegin(); node != faces.NodesEnd(); ++node)
	  {
		  const Vector& local_parameters = node->GetValue(LOCAL_PARAMETERS);
		  file << separator << "[" << node->Id() << ", " <<  node->GetValue(FACE_BREP_ID) << ",[" << local_parameters[0] << ", " << local_parameters[1] << "]]";
		  separator = ",";
	  }
	  file << "]}";
	  file.close();
  }


  double BrepModelGeometryReader::ReadModelTolerance()
  {
	  return m_cad_geometry_in_json["tolerances"]["model_tolerance"].GetDouble();
  }
  //Parameters BrepModelGeometryReader::WriteGeometry(ModelPart& model_part)
  //{
	 // std::cout << "\n> Start write CAD geometry..." << std::endl;
	 // Paramaters geometry;
	 // geometry.
  //}
  std::vector<BrepModel> BrepModelGeometryReader::ReadGeometry(ModelPart& model_part)
  {
      std::cout << "\n> Start reading CAD geometry" << std::endl;


      std::vector<BrepModel> r_brep_model_vector;

      for (int brep_i = 0; brep_i < m_cad_geometry_in_json["breps"].size(); brep_i++)
      {
          Parameters brep_json = m_cad_geometry_in_json["breps"][brep_i];
          unsigned int brep_brep_id = brep_json["brep_id"].GetInt();

          std::vector<BrepFace>   faces_vector;
          std::vector<BrepEdge>   edges_vector;
          std::vector<BrepVertex> vertices_vector;

          // loop over faces
          for (int i = 0; i < brep_json["faces"].size(); i++)
          {

              /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              // 1. Step: faces
              /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              unsigned int face_id = brep_json["faces"][i]["brep_id"].GetInt();

              //model_part.CreateSubModelPart("FACE_" + std::to_string(face_id));
              ModelPart& sub_model_part_face = model_part.CreateSubModelPart("FACE_" + std::to_string(face_id));

              //sub_model_part_face.CreateSubModelPart("FACE_" + std::to_string(face_id) + "_CPS");
              ModelPart& sub_model_part_face_cp = sub_model_part_face.CreateSubModelPart("FACE_" + std::to_string(face_id) + "_CPS");

              //std::cout << "> Reading face " << face_id << "..." << std::endl;

              bool is_trimmed = brep_json["faces"][i]["surface"]["is_trimmed"].GetBool();
              bool is_rational = brep_json["faces"][i]["surface"]["is_rational"].GetBool();

              // Variables needed later
              unsigned int length_u_vector = brep_json["faces"][i]["surface"]["knot_vectors"][0].size();
              unsigned int length_v_vector = brep_json["faces"][i]["surface"]["knot_vectors"][1].size();
              Vector knot_vector_u = ZeroVector(length_u_vector + 2);
              Vector knot_vector_v = ZeroVector(length_v_vector + 2);
              unsigned int p;
              unsigned int q;
              IntVector control_points_ids;

              // read and store knot_vector_u
              knot_vector_u(0) = brep_json["faces"][i]["surface"]["knot_vectors"][0][0].GetDouble();
              for (int u_idx = 0; u_idx < length_u_vector; u_idx++)
              {
                  knot_vector_u(u_idx + 1) = brep_json["faces"][i]["surface"]["knot_vectors"][0][u_idx].GetDouble();
                  //knot_vector_u.insert_element(-1, knot);
              }
              knot_vector_u[length_u_vector + 1] = brep_json["faces"][i]["surface"]["knot_vectors"][0][length_u_vector - 1].GetDouble();

              // read and store knot_vector_v
              knot_vector_v(0) = brep_json["faces"][i]["surface"]["knot_vectors"][1][0].GetDouble();
              for (int v_idx = 0; v_idx < length_v_vector; v_idx++)
              {
                  knot_vector_v(v_idx + 1) = brep_json["faces"][i]["surface"]["knot_vectors"][1][v_idx].GetDouble();
                  //knot_vector_v.insert_element(-1, knot);
              }
              knot_vector_v[length_v_vector + 1] = brep_json["faces"][i]["surface"]["knot_vectors"][1][length_v_vector - 1].GetDouble();

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
                  sub_model_part_face_cp.pGetNode(cp_id)->SetValue(CONTROL_POINT_WEIGHT, w);
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

                  // Loop over all curves
                  for (int edge_idx = 0; edge_idx < boundary_dict[loop_idx]["trimming_curves"].size(); edge_idx++)
                  {
                      unsigned int curve_index = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["trim_index"].GetInt();
                      //loop.push_back(curve_index);

                      bool curve_direction = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["curve_direction"].GetBool();
                      // Variables needed later
                      unsigned int length_u_vector = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["knot_vector"].size();
                      Vector boundary_knot_vector_u = ZeroVector(length_u_vector + 2);
                      unsigned int boundary_p;
                      std::vector<array_1d<double, 4>> boundary_control_points;
                      Vector active_range = ZeroVector(2);
                      // read and store knot_vector_u
                      boundary_knot_vector_u(0) = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["knot_vector"][0].GetDouble();
                      for (int u_idx = 0; u_idx < length_u_vector; u_idx++)
                      {
                          boundary_knot_vector_u(u_idx + 1) = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["knot_vector"][u_idx].GetDouble();
                      }
                      boundary_knot_vector_u[length_u_vector + 1] = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["knot_vector"][length_u_vector - 1].GetDouble();

                      // read and store polynamial degree p and q
                      boundary_p = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["degree"].GetInt();

                      // read and store control_points
                      for (int cp_idx = 0; cp_idx < boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"].size(); cp_idx++)
                      {
                          //unsigned int cp_id = extractInt(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][0]);
                          array_1d<double, 4> control_point;
                          control_point[0] = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][0].GetDouble();
                          control_point[1] = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1].GetDouble();
                          control_point[2] = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][2].GetDouble();
                          control_point[3] = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][3].GetDouble();

                          boundary_control_points.push_back(control_point);
                      }
                      active_range(0) = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["active_range"][0].GetDouble();
                      active_range(1) = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["active_range"][1].GetDouble();

                      // Create and store edge
                      BrepTrimmingCurve new_boundary_curve(curve_index, curve_direction, boundary_knot_vector_u, boundary_p, boundary_control_points, active_range);

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
                      unsigned int curve_index = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["trim_index"].GetInt();
                      //loop.push_back(curve_index);

                      bool curve_direction = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["curve_direction"].GetBool();
                      // Variables needed later
                      unsigned int length_u_vector = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["knot_vector"].size();
                      Vector boundary_knot_vector_u = ZeroVector(length_u_vector + 2);
                      unsigned int boundary_p;
                      std::vector<array_1d<double, 4>> boundary_control_points;
                      Vector active_range = ZeroVector(2);
                      // read and store knot_vector_u
                      boundary_knot_vector_u(0) = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["knot_vector"][0].GetDouble();
                      for (int u_idx = 0; u_idx < length_u_vector; u_idx++)
                      {
                          boundary_knot_vector_u(u_idx + 1) = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["knot_vector"][u_idx].GetDouble();
                      }
                      boundary_knot_vector_u[length_u_vector + 1] = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["knot_vector"][length_u_vector - 1].GetDouble();

                      // read and store polynamial degree p and q
                      boundary_p = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["degree"].GetInt();

                      // read and store control_points
                      for (int cp_idx = 0; cp_idx < embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"].size(); cp_idx++)
                      {
                          //unsigned int cp_id = extractInt(embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][0]);
                          array_1d<double, 4> control_point;
                          control_point[0] = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][0].GetDouble();
                          control_point[1] = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1].GetDouble();
                          control_point[2] = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][2].GetDouble();
                          control_point[3] = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][3].GetDouble();

                          boundary_control_points.push_back(control_point);
                      }
                      active_range(0) = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["active_range"][0].GetDouble();
                      active_range(1) = embedded_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["active_range"][1].GetDouble();

                      // Create and store edge
                      BrepTrimmingCurve new_boundary_curve(curve_index, curve_direction, boundary_knot_vector_u, boundary_p, boundary_control_points, active_range);

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
              Kratos::shared_ptr<ModelPart> sub_model_part_face_ptr = std::shared_ptr<ModelPart>(sub_model_part_face.pGetSubModelPart("FACE_" + std::to_string(face_id) + "_CPS"));
              // create face
              BrepFace face(face_id, is_trimmed, is_rational, trimming_loops, embedded_loops, embedded_points, knot_vector_u, knot_vector_v, p, q, 
                  control_points_ids, sub_model_part_face_ptr);
              faces_vector.push_back(face);
          }

          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          // Edges
          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          for (int i = 0; i < brep_json["edges"].size(); i++)
          {
              // For better reading
              Parameters edge_dict(brep_json["edges"][i]);

              std::vector<BrepEdge::Topology> brep_edge_topology_vector;

              unsigned int edge_id = edge_dict["brep_id"].GetInt();

              //model_part.CreateSubModelPart("FACE_" + std::to_string(face_id));
              ModelPart& sub_model_part_edge = model_part.CreateSubModelPart("EDGE_" + std::to_string(edge_id));

              std::cout << "> Reading edge " << edge_id << "..." << std::endl;

              /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              // 3d curve
              /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

              unsigned int degree = edge_dict["3d_curve"]["degree"].GetInt();

              // read and store active_range
              Vector active_range(2);
              active_range(0) = edge_dict["3d_curve"]["active_range"][0].GetDouble();
              active_range(1) = edge_dict["3d_curve"]["active_range"][1].GetDouble();

              // read and store knot_vector
              double length_knot_vector = edge_dict["3d_curve"]["knot_vector"].size();
              Vector knot_vector = ZeroVector(length_knot_vector + 2);
              knot_vector(0) = edge_dict["3d_curve"]["knot_vector"][0].GetDouble();
              for (int u_idx = 0; u_idx < length_knot_vector; u_idx++)
              {
                  knot_vector(u_idx + 1) = edge_dict["3d_curve"]["knot_vector"][u_idx].GetDouble();
              }
              knot_vector(length_knot_vector + 1) = edge_dict["3d_curve"]["knot_vector"][length_knot_vector - 1].GetDouble();


              IntVector control_points_ids_edge;
              // read and store control_points
              // Control points in each patch get a global as well as a mapping matrix id
              // brep Id: Unique Id for each control point (given by json-file)
              // mapping matrix id: specifies position in global mapping matrix (given by numbering of control points)
              for (int cp_idx = 0; cp_idx < edge_dict["3d_curve"]["control_points"].size(); cp_idx++)
              {
                  unsigned int cp_id = edge_dict["3d_curve"]["control_points"][cp_idx][0].GetInt();
                  double x = edge_dict["3d_curve"]["control_points"][cp_idx][1][0].GetDouble();
                  double y = edge_dict["3d_curve"]["control_points"][cp_idx][1][1].GetDouble();
                  double z = edge_dict["3d_curve"]["control_points"][cp_idx][1][2].GetDouble();
                  double w = edge_dict["3d_curve"]["control_points"][cp_idx][1][3].GetDouble();

                  control_points_ids_edge.push_back(cp_id);
                  //sub_model_part_edge->CreateNewNode(cp_id, x, y, z);
                  //sub_model_part_edge->GetNode(cp_id).SetValue(CONTROL_POINT_WEIGHT, w);
              }

              /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              // trimming ranges
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
              // topology
              /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              for (int j = 0; j < edge_dict["topology"].size(); j++)
              {
                  unsigned int face_id = edge_dict["topology"][j]["face_id"].GetInt();
                  unsigned int trim_index = edge_dict["topology"][j]["trim_index"].GetInt();
                  bool relative_direction = edge_dict["topology"][j]["relative_direction"].GetBool();

                  BrepEdge::Topology trim(face_id, trim_index, relative_direction);
                  brep_edge_topology_vector.push_back(trim);
              }

              std::cout << "> Finished reading edge " << edge_id << "..." << std::endl;
              Kratos::shared_ptr<ModelPart> sub_model_part_edge_ptr = std::shared_ptr<ModelPart>(model_part.pGetSubModelPart("EDGE_" + std::to_string(edge_id)));
              BrepEdge edge(edge_id, brep_edge_topology_vector, brep_trimming_range_vector,
                  degree, knot_vector, active_range, control_points_ids_edge, sub_model_part_edge_ptr);

              std::cout << "> Finished creating edge " << edge_id << "..." << std::endl;
              edges_vector.push_back(edge);
          }
          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          // Vertices
          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

          if (brep_json.Has("vertices"))
          {
              for (int i = 0; i < brep_json["vertices"].size(); i++)
              {
                  // For better reading
                  Parameters vertex_dict(brep_json["vertices"][i]);

                  unsigned int vertex_id = vertex_dict["brep_id"].GetInt();

                  std::cout << "> Start reading vertex " << vertex_id << "..." << std::endl;

                  int cp_id = 0;// vertex_dict["coordinates"][0].GetInt();
                  Vector coords = ZeroVector(4);
                  //coords(0) = vertex_dict["coordinates"][1][0].GetDouble();
                  //coords(1) = vertex_dict["coordinates"][1][1].GetDouble();
                  //coords(2) = vertex_dict["coordinates"][1][2].GetDouble();
                  //coords(3) = vertex_dict["coordinates"][1][3].GetDouble();

                  std::vector<BrepVertex::Topology> brep_vertices_topology_vector;
                  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  // topology
                  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  for (int j = 0; j < vertex_dict["topology"].size(); j++)
                  {
                      unsigned int face_id = vertex_dict["topology"][j]["brep_id"].GetInt();
                      unsigned int trim_index = vertex_dict["topology"][j]["trim_index"].GetInt();

                      BrepVertex::Topology trim(face_id, trim_index);
                      brep_vertices_topology_vector.push_back(trim);
                  }
                  std::cout << "> Finished reading vertex " << vertex_id << "..." << std::endl;

                  BrepVertex vertex(vertex_id, brep_vertices_topology_vector, cp_id, coords);
                  vertices_vector.push_back(vertex);
                  std::cout << "> Finished creating vertex " << vertex_id << "..." << std::endl;
              }
          }
          std::cout << "Before creating BrepModel" << std::endl;

          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          // 3. Step: Create BrepModel
          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          BrepModel brep(brep_brep_id, faces_vector, edges_vector, vertices_vector);

          r_brep_model_vector.push_back(brep);// [brep_i] = &brep;
      }
      std::cout << "\n> Finished reading CAD geometry" << std::endl;
      return r_brep_model_vector;
  }

// --------------------------------------------------------------------------

//Constructor
BrepModelGeometryReader::BrepModelGeometryReader(Parameters& cad_geometry_in_json)
  : m_cad_geometry_in_json(cad_geometry_in_json)
{
}
//Destructor
BrepModelGeometryReader::~BrepModelGeometryReader()
{}

}  // namespace Kratos.
