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
#include "iga_application_variables.h"

namespace Kratos
{

///@name Stages
///@{

void ImportNurbsSbmModeler::SetupGeometryModel(){
    // Get bounding box physical space.

    KRATOS_ERROR_IF_NOT(mParameters.Has("input_filename"))
        << "::[ImportNurbsSbmModeler]:: Missing \"input_filename\" section." << std::endl;

    KRATOS_ERROR_IF_NOT(mParameters.Has("model_part_name"))
        << "::[ImportNurbsSbmModeler]:: Missing \"model_part_name\" section." << std::endl;
    
    KRATOS_ERROR_IF_NOT(mParameters.Has("link_layer_to_condition_name"))
        << "::[ImportNurbsSbmModeler]:: Missing \"link_layer_to_condition_name\" section." << std::endl;

    // Create the NURBS skin geometries from the json file
    const std::string input_filename = mParameters["input_filename"].GetString();
    const std::string model_part_name = mParameters["model_part_name"].GetString();
    auto layer_condition_name = mParameters["link_layer_to_condition_name"];

    ModelPart& skin_model_part_initial = mpModel->CreateModelPart(model_part_name);
    skin_model_part_initial.CreateNewProperties(0);
    const Parameters nurbs_skin_parameters = ReadParamatersFile(input_filename);
    
    // Get the number of boundary NURBS
    SizeType n_boundary_curves = nurbs_skin_parameters["Lines"].size();

    SizeType last_node_id = 0;
    for (IndexType i_curve = 0; i_curve < n_boundary_curves; i_curve++)
    {
        // Check NURBS geometry features
        CheckNurbsFeatures(nurbs_skin_parameters["Lines"][i_curve]);
 
        SizeType n_control_points = nurbs_skin_parameters["Lines"][i_curve]["Weights"].size();
        auto control_point_coordinates = nurbs_skin_parameters["Lines"][i_curve]["CPCoordinates"];
        int polynomial_degree = nurbs_skin_parameters["Lines"][i_curve]["pDegree"].GetInt();
        Vector knot_vector = nurbs_skin_parameters["Lines"][i_curve]["knotVector"].GetVector();
        Vector weights_vector = nurbs_skin_parameters["Lines"][i_curve]["Weights"].GetVector();
        std::string layer_name = nurbs_skin_parameters["Lines"][i_curve]["Layer"].GetString();
        std::string condition_name;

        // find associated condition_name
        bool layer_exits = false;
        for (IndexType i_layer = 0; i_layer < layer_condition_name.size(); i_layer++){
            if (layer_condition_name[i_layer]["layer_name"].GetString() == layer_name)
            {
                condition_name = layer_condition_name[i_layer]["condition_name"].GetString();
                layer_exits = true;
                break;
            }
        }
        KRATOS_ERROR_IF_NOT(layer_exits)
            << "::[ImportNurbsSbmModeler]:: geometry layer \"" << layer_name << "\" does not match any layer_condition_name" << std::endl
            << "layer_condition_name availables are: " << std::endl << layer_condition_name << std::endl;
        
        // store the control points
        PointerVector<Node> control_points;
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
        
        // link the boundary condition and layer name to the nurbs curve 
        p_curve->SetValue(CONDITION_NAME, condition_name);
        p_curve->SetValue(IDENTIFIER, layer_name);

        // add nurbs curve to the model part 
        p_curve->SetId(i_curve);
        skin_model_part_initial.AddGeometry(p_curve);
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

void ImportNurbsSbmModeler::CheckNurbsFeatures(const Parameters& rNurbsCurveParameters) {
    KRATOS_ERROR_IF_NOT(rNurbsCurveParameters.Has("CPCoordinates")) << "Control Points not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsCurveParameters.Has("pDegree")) << "Curve degree not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsCurveParameters.Has("knotVector")) << "Curve Knot Vector not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsCurveParameters.Has("Weights")) << "Curve weights not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsCurveParameters.Has("Layer")) << "Curve layer not defined in the NURBS geometry file." << std::endl;
}
    
} // end namespace kratos
