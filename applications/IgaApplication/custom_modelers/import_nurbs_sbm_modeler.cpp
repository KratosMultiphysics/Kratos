//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi
//

// Project includes
#include "includes/define.h"
#include "import_nurbs_sbm_modeler.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities.h"
#include "iga_application_variables.h"

namespace Kratos
{

    ///@name Stages
    ///@{

    void ImportNurbsSbmModeler::SetupGeometryModel(){
        // Get bounding box physical space.

        KRATOS_ERROR_IF_NOT(mParameters.Has("input_filename"))
            << "ImportNurbsSbmModeler: Missing \"input_filename\" section." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("model_part_name"))
            << "ImportNurbsSbmModeler: Missing \"model_part_name\" section." << std::endl;


        // Create the NURBS skin geometries from the json file
        const std::string input_filename = mParameters["input_filename"].GetString();
        const std::string model_part_name = mParameters["model_part_name"].GetString();
        auto layer_condition_name = mParameters["associate_layer_with_condition_name"];

        ModelPart& skin_model_part_initial = mpModel->CreateModelPart(model_part_name);
        skin_model_part_initial.CreateNewProperties(0);

        const Parameters nurbs_skin_parameters = ReadParamatersFile(input_filename);
        
        SizeType n_boundary_curves = nurbs_skin_parameters["Lines"].size();

        // maybe for the future to check the orientation of the boundary lines
        // std::vector<NurbsCurveGeometryPointerType> boundary_curves(n_boundary_curves);

        SizeType last_node_id = 0;
        for (IndexType i_curve = 0; i_curve < n_boundary_curves; i_curve++)
        {
            SizeType n_control_points = nurbs_skin_parameters["Lines"][i_curve]["Weights"].size();
            auto control_point_coordinates = nurbs_skin_parameters["Lines"][i_curve]["CPCoordinates"];

            PointerVector<Node> control_points;
            int polynomial_degree = nurbs_skin_parameters["Lines"][i_curve]["pDegree"].GetInt();
            Vector knot_vector = nurbs_skin_parameters["Lines"][i_curve]["knotVector"].GetVector();
            Vector weights_vector = nurbs_skin_parameters["Lines"][i_curve]["Weights"].GetVector();

            std::string layer_name = nurbs_skin_parameters["Lines"][i_curve]["Layer"].GetString();
            std::string condition_name;

            // find associated condition_name
            for (IndexType i_layer = 0; i_layer < layer_condition_name.size(); i_layer++){
                if (layer_condition_name[i_layer]["layer_name"].GetString() == layer_name)
                {
                    condition_name = layer_condition_name[i_layer]["condition_name"].GetString();
                }
            } 

            for (IndexType i_cp = 0; i_cp < n_control_points; i_cp++)
            {
                Vector current_cp_coordinates = control_point_coordinates[i_cp].GetVector();
                control_points.push_back(Kratos::make_intrusive<NodeType>(0, current_cp_coordinates));
            }

            NurbsCurveGeometryPointerType p_curve(new NurbsCurveGeometry<2, PointerVector<Node>>(
                                                            control_points,
                                                            polynomial_degree,
                                                            knot_vector, 
                                                            weights_vector)); 

            p_curve->SetValue(CONDITION_NAME, condition_name);
            p_curve->SetValue(IDENTIFIER, layer_name);

            // add curve to the model part 
            p_curve->SetId(i_curve);
            skin_model_part_initial.AddGeometry(p_curve);

            // // first point
            // CoordinatesArrayType point_coords(3);
            // Vector local_coord = ZeroVector(3);
            // p_curve->GlobalCoordinates(point_coords, local_coord);

            // // check the first point of the curve
            // if (i_curve == 0) 
            // {
            //     last_node_id++;
            //     Node::Pointer node = new Node(last_node_id, point_coords[0], point_coords[1], point_coords[2]);
            //     skin_model_part_initial.CreateNewNode(last_node_id, point_coords[0], point_coords[1], point_coords[2]); 
            // } else 
            // {
            //     KRATOS_ERROR_IF(norm_2(point_coords - skin_model_part_initial.GetNode(last_node_id)) > 1e-15)
            //                     << "ImportNurbsSbmModeler: error in the json NURBS file. The boundary curves" 
            //                     << " are not correctly ordered." << std::endl;;
            // }

            // const SizeType first_node_id = last_node_id;
            // // add the specified number of points
            // for (int i = 1; i < n_initial_points_for_side-1; i++)
            // {
            //     local_coord[0] = (double) i/(n_initial_points_for_side-1);
            //     p_curve->GlobalCoordinates(point_coords, local_coord);
            //     last_node_id++;
            //     // Node::Pointer node = new Node(last_node_id, point_coords[0], point_coords[1], point_coords[2]);

            //     // // compute normal and store value to the point
            //     // std::vector<CoordinatesArrayType> global_space_derivatives;
            //     // SizeType derivative_order = 1;
            //     // p_curve->GlobalSpaceDerivatives(global_space_derivatives, local_coord, derivative_order);
            //     // CoordinatesArrayType tangent_vector = global_space_derivatives[1];
            //     // double tangent_magnitude = norm_2(tangent_vector);
            //     // tangent_vector /= tangent_magnitude;
            //     // Vector normal_vector = ZeroVector(3);
            //     // normal_vector[0] = tangent_vector[1];
            //     // normal_vector[1] = -tangent_vector[0];

            //     // // node->SetValue(NORMAL, node)
            //     // skin_model_part_initial.AddNode(node);

            //     //-------------------------------------------------
            //     skin_model_part_initial.CreateNewNode(last_node_id, point_coords[0], point_coords[1], point_coords[2]); 
            //     skin_model_part_initial.CreateNewCondition("LineCondition2D2N", last_node_id-1, {{last_node_id-1, last_node_id}}, 0);
            // }
            // // check the last point of the curve
            // local_coord[0] = 1.0;
            // p_curve->GlobalCoordinates(point_coords, local_coord);
            // if (i_curve == n_boundary_curves-1) 
            // {
            //     skin_model_part_initial.CreateNewCondition("LineCondition2D2N", last_node_id, {{last_node_id, 1}}, 0);
            //     KRATOS_ERROR_IF(norm_2(point_coords - skin_model_part_initial.GetNode(1)) > 1e-15)
            //                     << "ImportNurbsSbmModeler: error in the json NURBS file. The skin boundary is not closed!" << std::endl;;
            // } else 
            // {
            //     last_node_id++;
            //     skin_model_part_initial.CreateNewNode(last_node_id, point_coords[0], point_coords[1], point_coords[2]);
            //     skin_model_part_initial.CreateNewCondition("LineCondition2D2N", last_node_id-1, {{last_node_id-1, last_node_id}}, 0);
            // }    

            // // create fictitious element to keep track of which conditions are taken from which curve
            // std::vector<ModelPart::IndexType> elem_nodes{first_node_id, last_node_id};
            // skin_model_part_initial.CreateNewElement("Element2D2N", skin_model_part_initial.GetRootModelPart().Elements().size()+1, elem_nodes, 0);
        } 

    }
    ///@}



    /// Reads in a json formatted file and returns its KratosParameters instance.
    Parameters ImportNurbsSbmModeler::ReadParamatersFile(
        const std::string& rDataFileName) const
    {
        // Check if rDataFileName ends with ".cad.json" and add it if needed.
        const std::string data_file_name = (rDataFileName.compare(rDataFileName.size() - 5, 5, ".json") != 0)
            ? rDataFileName + ".json"
            : rDataFileName;

        std::ifstream infile(data_file_name);
        KRATOS_ERROR_IF_NOT(infile.good()) << "Nurbs geometry file: "
            << data_file_name << " cannot be found." << std::endl;
        KRATOS_INFO_IF("ReadParamatersFile", mEchoLevel > 3)
            << "Reading file: \"" << data_file_name << "\"" << std::endl;

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    }


    /* Reads, and returns a Pointer to Node.
    * Input needs to be a Parameter object:
    * [id, [x, y, z, weight]]
    */
    // static Node::Pointer ReadAndCreateNode(
    //     const Parameters rParameters,
    //     ModelPart& rModelPart,
    //     SizeType EchoLevel = 0)
    // {
    //     SizeType number_of_entries = rParameters.size();
    //     KRATOS_ERROR_IF((number_of_entries != 2))
    //         << "Control points as Node need to be provided in following structure: [id, [x, y, z, weight]]"
    //         << std::endl;

    //     IndexType id = rParameters[0].GetInt();
    //     Vector cp = rParameters[1].GetVector();

    //     return rModelPart.CreateNewNode(id, cp[0], cp[1], cp[2]);
    // }
} // end namespace kratos
