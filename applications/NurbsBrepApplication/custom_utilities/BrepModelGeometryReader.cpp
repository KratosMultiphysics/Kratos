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
  std::vector<BrepModel> BrepModelGeometryReader::ReadGeometry(ModelPart& model_part)
  {
    std::cout << "\n> Start reading CAD geometry..." << std::endl;

    std::vector<BrepModel> r_brep_model_vector;

    for (int brep_i = 0; brep_i < m_cad_geometry_in_json["breps"].size(); brep_i++)
    {
      Parameters& brep_json = m_cad_geometry_in_json["breps"][brep_i];
      unsigned int brep_brep_id = brep_json["brep_id"].GetInt();

      BrepFacesVector faces_vector;
      BrepEdgesVector edges_vector;

      // loop over faces
      for (int i = 0; i < brep_json["faces"].size(); i++)
      {

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 1. Step: faces
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        unsigned int face_id = brep_json["faces"][i]["brep_id"].GetInt();

        ModelPart& sub_model_part_face = model_part.CreateSubModelPart("FACE_" + std::to_string(face_id));
        ModelPart& sub_model_part_face_cp = sub_model_part_face.CreateSubModelPart("FACE_" + std::to_string(face_id) + "_CPS");

        std::cout << "> Reading face " << face_id << "..." << std::endl;

        // Variables needed later
        unsigned int length_u_vector = brep_json["faces"][i]["surface"]["knot_vectors"][0].size();
        unsigned int length_v_vector = brep_json["faces"][i]["surface"]["knot_vectors"][1].size();
        Vector knot_vector_u = ZeroVector(length_u_vector+2);
        Vector knot_vector_v = ZeroVector(length_v_vector+2);
        unsigned int p;
        unsigned int q;
        IntVector control_points_ids;

        // read and store knot_vector_u
        knot_vector_u(0) = brep_json["faces"][i]["surface"]["knot_vectors"][0][0].GetDouble();
        for (int u_idx = 0; u_idx < length_u_vector; u_idx++)
        {
          knot_vector_u(u_idx+1) = brep_json["faces"][i]["surface"]["knot_vectors"][0][u_idx].GetDouble();
          //knot_vector_u.insert_element(-1, knot);
        }
        knot_vector_u[length_u_vector+1] = brep_json["faces"][i]["surface"]["knot_vectors"][0][length_u_vector-1].GetDouble();

        // read and store knot_vector_v
        knot_vector_v(0) = brep_json["faces"][i]["surface"]["knot_vectors"][1][0].GetDouble();
        for (int v_idx = 0; v_idx < length_v_vector; v_idx++)
        {
          knot_vector_v(v_idx+1) = brep_json["faces"][i]["surface"]["knot_vectors"][1][v_idx].GetDouble();
          //knot_vector_v.insert_element(-1, knot);
        }
        knot_vector_v[length_v_vector + 1] = brep_json["faces"][i]["surface"]["knot_vectors"][1][length_v_vector - 1].GetDouble();

        // read and store polynamial degree p and q
        p = brep_json["faces"][i]["surface"]["degrees"][0].GetInt();
        q = brep_json["faces"][i]["surface"]["degrees"][1].GetInt();


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
          sub_model_part_face_cp.GetNode(cp_id).SetValue(CONTROL_POINT_WEIGHT, w);
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 2. step: boundary loops
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        BrepTrimmingLoopVector trimming_loops;
        //TrimmingCurveVector trimming_curves;

        // For better reading
        Parameters boundary_dict(brep_json["faces"][i]["boundary_loops"]);

        // Loop over all boundary loops
        for (int loop_idx = 0; loop_idx < boundary_dict.size(); loop_idx++)
        {
          BrepTrimmingCurveVector loop_curves;

          // Loop over all curves
          for (int edge_idx = 0; edge_idx < boundary_dict[loop_idx]["trimming_curves"].size(); edge_idx++)
          {
            unsigned int curve_index = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["trim_index"].GetInt();
            //loop.push_back(curve_index);

            bool curve_direction = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["curve_direction"].GetBool();
            // Variables needed later
            unsigned int length_u_vector = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["u_vec"].size();
            Vector boundary_knot_vector_u = ZeroVector(length_u_vector+2);
            unsigned int boundary_p;
            std::vector<array_1d<double, 4>> boundary_control_points;
            Vector active_range = ZeroVector(2);
            // read and store knot_vector_u
            boundary_knot_vector_u(0) = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["u_vec"][0].GetDouble();
            for (int u_idx = 0; u_idx < length_u_vector; u_idx++)
            {
              boundary_knot_vector_u(u_idx+1) = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["u_vec"][u_idx].GetDouble();
            }
            boundary_knot_vector_u[length_u_vector + 1] = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["u_vec"][length_u_vector - 1].GetDouble();

            // read and store polynamial degree p and q
            boundary_p = boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["degrees"].GetInt();

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

        // create face
        BrepFace face(face_id, trimming_loops, knot_vector_u, knot_vector_v, p, q, control_points_ids, sub_model_part_face);
        faces_vector.push_back(face);
      }

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Edges
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for (int i = 0; i < brep_json["edges"].size(); i++)
      {
        BrepFaceTrimVector face_trims_vector;

        unsigned int edge_id = brep_json["edges"][i]["brep_id"].GetInt();
        Vector first_boundary_vertex = ZeroVector(3);
        Vector second_boundary_vertex = ZeroVector(3);
        ParameterVector boundary_vertices;
        for (int j = 0; j < 3; j++)
        {
          first_boundary_vertex(j) = brep_json["edges"][i]["boundary_vertices"][0][j].GetDouble();
          second_boundary_vertex(j) = brep_json["edges"][i]["boundary_vertices"][1][j].GetDouble();
        }
        boundary_vertices.push_back(first_boundary_vertex);
        boundary_vertices.push_back(second_boundary_vertex);
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // FaceTrims
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        for (int j = 0; j < brep_json["edges"][i]["face_trims"].size(); j++)
        {
          unsigned int face_id = brep_json["edges"][i]["face_trims"][j]["face_id"].GetInt();
          unsigned int trim_index = brep_json["edges"][i]["face_trims"][j]["trim_index"].GetInt();

          bool relative_direction = brep_json["edges"][i]["face_trims"][j]["relative_direction"].GetBool();

          BrepFaceTrim trim(face_id, trim_index, relative_direction);
          face_trims_vector.push_back(trim);
        }
        BrepEdge edge(edge_id, boundary_vertices, face_trims_vector);
        edges_vector.push_back(edge);
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 3. Step: Create BrepModel
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      BrepModel brep(brep_brep_id, faces_vector, edges_vector);
      
      r_brep_model_vector.push_back(brep);// [brep_i] = &brep;
    }
    std::cout << "\n> Finished reading CAD geometry..." << std::endl;
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


