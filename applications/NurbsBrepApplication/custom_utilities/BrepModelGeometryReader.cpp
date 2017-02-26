//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
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
  void BrepModelGeometryReader::ReadGeometry(BrepModelVector& r_brep_model_vector, ModelPart& model_part)
  {
    std::cout << "\n> Start reading CAD geometry..." << std::endl;

    FacesVector faces_vector;
    EdgesVector edges_vector;

    // loop over faces
    for (int i = 0; i < boost::python::len(m_cad_geometry_in_json["faces"]); i++)
    {
      
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 1. Step: faces
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      unsigned int face_id = extractInt(m_cad_geometry_in_json["faces"][i]["brep_id"]);

      std::string sub_model_part_face_name = "FACE_" + std::to_string(face_id);
      std::string sub_model_part_name = "FACE_" + std::to_string(face_id) + "_CPS";
      //std::string sub_model_part_name_internal = "FACE_" + std::to_string(face_id) + "_INTERNAL_CPS";

      std::cout << "Sub_model_part_name: " << sub_model_part_name << std::endl;
      //std::cout << "Sub_model_part_name_internal: " << sub_model_part_name_internal << std::endl;

      ModelPart& sub_model_part_face = model_part.CreateSubModelPart(sub_model_part_face_name);
      ModelPart& sub_model_part_face_cp = sub_model_part_face.CreateSubModelPart(sub_model_part_name);
      //ModelPart& sub_model_part_face_internal = sub_model_part_face.CreateSubModelPart(sub_model_part_name_internal);
      
      std::cout << "> Reading face " << face_id << "..." << std::endl;

      // Variables needed later
      DoubleVector knot_vector_u;
      DoubleVector knot_vector_v;
      int p;
      int q;
      IntVector control_points_ids;

      // read and store knot_vector_u
      for (int u_idx = 0; u_idx < boost::python::len(m_cad_geometry_in_json["faces"][i]["surface"][0]["knot_vectors"][0]); u_idx++)
      {
        double knot = extractDouble(m_cad_geometry_in_json["faces"][i]["surface"][0]["knot_vectors"][0][u_idx]);
        knot_vector_u.push_back(knot);
      }

      // read and store knot_vector_v
      for (int v_idx = 0; v_idx < boost::python::len(m_cad_geometry_in_json["faces"][i]["surface"][0]["knot_vectors"][1]); v_idx++)
      {
        double knot = extractDouble(m_cad_geometry_in_json["faces"][i]["surface"][0]["knot_vectors"][1][v_idx]);
        knot_vector_v.push_back(knot);
      }

      // read and store polynamial degree p and q
      p = extractInt(m_cad_geometry_in_json["faces"][i]["surface"][0]["degrees"][0]);
      q = extractInt(m_cad_geometry_in_json["faces"][i]["surface"][0]["degrees"][1]);

    // read and store control_points
      // Control points in each patch get a global as well as a mapping matrix id
      // brep Id: Unique Id for each control point (given by json-file)
      // mapping matrix id: specifies position in global mapping matrix (given by numbering of control points)
      for (int cp_idx = 0; cp_idx < boost::python::len(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"]); cp_idx++)
      {
        unsigned int cp_id = extractInt(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][0]);
        double x = extractDouble(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][0]);
        double y = extractDouble(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][1]);
        double z = extractDouble(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][2]);
        double w = extractDouble(m_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][3]);

        control_points_ids.push_back(cp_id);

        sub_model_part_face_cp.CreateNewNode(cp_id, x, y, z);
        sub_model_part_face_cp.GetNode(cp_id).SetValue(CONTROL_POINT_WEIGHT, w);
        //ControlPoint new_control_point(x, y, z, w, cp_id);
        //model_part.CreateNewNode(cp_id, x, y, z);
        //control_points.push_back(new_control_point);
      }

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 2. step: boundary loops
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      TrimmingLoopVector trimming_loop;
      TrimmingCurveVector trimming_curves;

      // For better reading
      boost::python::list boundary_dict(m_cad_geometry_in_json["faces"][i]["boundary_loops"]);

      // Loop over all boundary loops
      for (int loop_idx = 0; loop_idx < boost::python::len(boundary_dict); loop_idx++)
      {
        IntVector loop;

        // Loop over all curves
        for (int edge_idx = 0; edge_idx < boost::python::len(boundary_dict[loop_idx]["trimming_curves"]); edge_idx++)
        {
          unsigned curve_index = extractInt(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["trim_index"]);

          loop.push_back(curve_index);

          // Variables needed later
          DoubleVector boundary_knot_vector_u;
          unsigned int boundary_p;
          std::vector<array_1d<double, 4>> boundary_control_points;
          DoubleVector boundary_vertices;

          // read and store knot_vector_u
          for (int u_idx = 0; u_idx < boost::python::len(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["u_vec"]); u_idx++)
          {
            double knot = extractDouble(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["u_vec"][u_idx]);
            boundary_knot_vector_u.push_back(knot);
          }

          // read and store polynamial degree p and q
          boundary_p = extractInt(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["degrees"]);

          // read and store control_points
          for (int cp_idx = 0; cp_idx < boost::python::len(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"]); cp_idx++)
          {
            //unsigned int cp_id = extractInt(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][0]);
            array_1d<double, 4> control_point;
            control_point[0] = extractDouble(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][0]);
            control_point[1] = extractDouble(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][1]);
            control_point[2] = extractDouble(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][2]);
            control_point[3] = extractDouble(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["parameter_curve"]["control_points"][cp_idx][3]);

            boundary_control_points.push_back(control_point);

            //sub_model_part_face_internal.CreateNewNode(cp_id, x, y, z);
            //sub_model_part_face_internal.GetNode(cp_id).SetValue(CONTROL_POINT_WEIGHT, w);
          }

          boundary_vertices.push_back(extractDouble(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["boundary_vertices"][0]));
          boundary_vertices.push_back(extractDouble(boundary_dict[loop_idx]["trimming_curves"][edge_idx]["boundary_vertices"][1]));

          //std::cout << "Test Point boundary vertices" << std::endl;

          // Create and store edge
          TrimmingCurve new_boundary_curve(curve_index, boundary_knot_vector_u, boundary_p, boundary_control_points, boundary_vertices);
          trimming_curves.push_back(new_boundary_curve);
        }

        // Read loop type
        //std::string loop_type = extractString(boundary_dict[loop_idx]["loop_type"]);
        //std::string Inner("Inner");
        //bool is_inner_loop = false;
        //if(loop_type.compare(Inner) == 0)
        //  is_inner_loop = true;

        // Create and store boundary loop
        trimming_loop.push_back(loop);
      }

      // create face
      Face face(face_id, trimming_curves, trimming_loop, knot_vector_u, knot_vector_v, p, q, control_points_ids);
      faces_vector.push_back(face);
    }

    // loop over edges
    for (int i = 0; i < boost::python::len(m_cad_geometry_in_json["edges"]); i++)
    {
      FaceTrimVector face_trims_vector;

      unsigned int edge_id = extractInt(m_cad_geometry_in_json["edges"][i]["brep_id"]);
      DoubleVector first_boundary_vertex;
      DoubleVector second_boundary_vertex;
      ParameterVector boundary_vertices;
      //std::cout << "Test Point boundary vertices 2" << std::endl;
      for (int j = 0; j < 3; j++)
      {
        first_boundary_vertex.push_back(extractDouble(m_cad_geometry_in_json["edges"][i]["boundary_vertices"][0][j]));

        //std::cout << "Test Point boundary vertices 2.1" << first_boundary_vertex[0] << std::endl;
        second_boundary_vertex.push_back(extractDouble(m_cad_geometry_in_json["edges"][i]["boundary_vertices"][1][j]));
      }
      boundary_vertices.push_back(first_boundary_vertex);
      boundary_vertices.push_back(second_boundary_vertex);

      for (int j = 0; j < boost::python::len(m_cad_geometry_in_json["edges"][i]["face_trims"]); j++)
      {
        unsigned int face_id = extractInt(m_cad_geometry_in_json["edges"][i]["face_trims"][j]["face_id"]);
        unsigned int trim_index = extractInt(m_cad_geometry_in_json["edges"][i]["face_trims"][j]["trim_index"]);

        //DoubleVector first_boundary_parameter;
        //DoubleVector second_boundary_parameter;
        //ParameterVector boundary_parameters;

        //std::cout << "Test Point boundary vertices 3" << std::endl;

        //for (int k = 0; k < 2; k++)
        //{
        //  first_boundary_parameter.push_back(extractDouble(m_cad_geometry_in_json["edges"][i]["face_trims"][j]["boundary_parameters"][0][k]));
        //  second_boundary_parameter.push_back(extractDouble(m_cad_geometry_in_json["edges"][i]["face_trims"][j]["boundary_parameters"][1][k]));
        //}
        //boundary_parameters.push_back(first_boundary_parameter);
        //boundary_parameters.push_back(second_boundary_parameter);

        bool relative_direction = extractBool(m_cad_geometry_in_json["edges"][i]["face_trims"][j]["relative_direction"]);
        
        FaceTrim trim(face_id, trim_index, relative_direction);
        face_trims_vector.push_back(trim);
      }
      Edge edge(edge_id, boundary_vertices, face_trims_vector);
      edges_vector.push_back(edge);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. Step: Create BrepModel
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    BrepModel brep(faces_vector, edges_vector);
    r_brep_model_vector.push_back(brep);

    std::cout << "\n> Finished reading CAD geometry..." << std::endl;
  }

// --------------------------------------------------------------------------

//Constructor
BrepModelGeometryReader::BrepModelGeometryReader(boost::python::dict cad_geometry_in_json)
  : m_cad_geometry_in_json(cad_geometry_in_json)
{
}
//Destructor
BrepModelGeometryReader::~BrepModelGeometryReader()
{}

}  // namespace Kratos.


