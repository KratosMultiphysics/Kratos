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
    const std::string geometry_type = mParameters.Has("geometry_type")
        ? mParameters["geometry_type"].GetString()
        : "Lines";

    KRATOS_ERROR_IF(geometry_type != "Lines" && geometry_type != "Surfaces")
        << "::[ImportNurbsSbmModeler]:: \"geometry_type\" must be either "
        << "\"Lines\" or \"Surfaces\", got \"" << geometry_type << "\"." << std::endl;

    ModelPart& skin_model_part_initial = mpModel->CreateModelPart(model_part_name);
    skin_model_part_initial.pGetProperties(0);
    const Parameters nurbs_skin_parameters = ReadParamatersFile(input_filename);

    KRATOS_WARNING_IF("::[ImportNurbsSbmModeler]::",
        !mParameters.Has("geometry_type") &&
        nurbs_skin_parameters.Has("Lines") &&
        nurbs_skin_parameters.Has("Surfaces") &&
        nurbs_skin_parameters["Lines"].size() > 0 &&
        nurbs_skin_parameters["Surfaces"].size() > 0)
        << "The input file contains both Lines and Surfaces, but \"geometry_type\" was not "
        << "specified. Importing Lines for backward compatibility. Set \"geometry_type\" to "
        << "\"Surfaces\" for a 3D NURBS skin." << std::endl;
    
    if (geometry_type == "Surfaces") {
        KRATOS_ERROR_IF_NOT(nurbs_skin_parameters.Has("Surfaces"))
            << "::[ImportNurbsSbmModeler]:: Missing \"Surfaces\" section in NURBS geometry file." << std::endl;

        const SizeType n_boundary_surfaces = nurbs_skin_parameters["Surfaces"].size();
        for (IndexType i_surface = 0; i_surface < n_boundary_surfaces; ++i_surface) {
            const Parameters surface_parameters = nurbs_skin_parameters["Surfaces"][i_surface];
            CheckNurbsSurfaceFeatures(surface_parameters);

            const Vector weights = surface_parameters["Weights"].GetVector();
            const auto control_point_coordinates = surface_parameters["CPCoordinates"];
            PointerVector<Node> control_points;
            for (IndexType i_cp = 0; i_cp < weights.size(); ++i_cp) {
                control_points.push_back(Kratos::make_intrusive<NodeType>(
                    0, control_point_coordinates[i_cp].GetVector()));
            }

            const std::string layer_name = surface_parameters["Layer"].GetString();
            auto p_surface = Kratos::make_shared<NurbsSurfaceGeometry<3, PointerVector<Node>>>(
                control_points,
                surface_parameters["pDegreeU"].GetInt(),
                surface_parameters["pDegreeV"].GetInt(),
                surface_parameters["knotVectors"]["Xi"].GetVector(),
                surface_parameters["knotVectors"]["Eta"].GetVector(),
                weights);

            p_surface->SetValue(CONDITION_NAME, GetConditionName(layer_name, layer_condition_name));
            p_surface->SetValue(IDENTIFIER, layer_name);
            p_surface->SetId(i_surface);
            skin_model_part_initial.AddGeometry(p_surface);
        }
        return;
    }

    KRATOS_ERROR_IF_NOT(nurbs_skin_parameters.Has("Lines"))
        << "::[ImportNurbsSbmModeler]:: Missing \"Lines\" section in NURBS geometry file." << std::endl;

    // Get the number of boundary NURBS
    SizeType n_boundary_curves = nurbs_skin_parameters["Lines"].size();

    for (IndexType i_curve = 0; i_curve < n_boundary_curves; i_curve++)
    {
        // Check NURBS geometry features
        CheckNurbsCurveFeatures(nurbs_skin_parameters["Lines"][i_curve]);
 
        SizeType n_control_points = nurbs_skin_parameters["Lines"][i_curve]["Weights"].size();
        auto control_point_coordinates = nurbs_skin_parameters["Lines"][i_curve]["CPCoordinates"];
        int polynomial_degree = nurbs_skin_parameters["Lines"][i_curve]["pDegree"].GetInt();
        Vector knot_vector = nurbs_skin_parameters["Lines"][i_curve]["knotVector"].GetVector();
        Vector weights_vector = nurbs_skin_parameters["Lines"][i_curve]["Weights"].GetVector();
        std::string layer_name = nurbs_skin_parameters["Lines"][i_curve]["Layer"].GetString();
        const std::string condition_name = GetConditionName(layer_name, layer_condition_name);
        
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

void ImportNurbsSbmModeler::CheckNurbsCurveFeatures(const Parameters& rNurbsCurveParameters) {
    KRATOS_ERROR_IF_NOT(rNurbsCurveParameters.Has("CPCoordinates")) << "Control Points not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsCurveParameters.Has("pDegree")) << "Curve degree not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsCurveParameters.Has("knotVector")) << "Curve Knot Vector not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsCurveParameters.Has("Weights")) << "Curve weights not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsCurveParameters.Has("Layer")) << "Curve layer not defined in the NURBS geometry file." << std::endl;
}

void ImportNurbsSbmModeler::CheckNurbsSurfaceFeatures(const Parameters& rNurbsSurfaceParameters) {
    KRATOS_ERROR_IF_NOT(rNurbsSurfaceParameters.Has("CPCoordinates")) << "Control Points not defined in the NURBS surface." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsSurfaceParameters.Has("pDegreeU")) << "Surface degree U not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsSurfaceParameters.Has("pDegreeV")) << "Surface degree V not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsSurfaceParameters.Has("knotVectors")) << "Surface knot vectors not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsSurfaceParameters["knotVectors"].Has("Xi")) << "Surface Xi knot vector not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsSurfaceParameters["knotVectors"].Has("Eta")) << "Surface Eta knot vector not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsSurfaceParameters.Has("Weights")) << "Surface weights not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF_NOT(rNurbsSurfaceParameters.Has("Layer")) << "Surface layer not defined in the NURBS geometry file." << std::endl;
    KRATOS_ERROR_IF(rNurbsSurfaceParameters["CPCoordinates"].size() != rNurbsSurfaceParameters["Weights"].size())
        << "Number of surface control points and weights do not match." << std::endl;
}

std::string ImportNurbsSbmModeler::GetConditionName(
    const std::string& rLayerName,
    const Parameters& rLayerConditionNames) const
{
    for (IndexType i_layer = 0; i_layer < rLayerConditionNames.size(); ++i_layer) {
        if (rLayerConditionNames[i_layer]["layer_name"].GetString() == rLayerName) {
            return rLayerConditionNames[i_layer]["condition_name"].GetString();
        }
    }
    KRATOS_ERROR << "::[ImportNurbsSbmModeler]:: geometry layer \"" << rLayerName
        << "\" does not match any layer_condition_name" << std::endl
        << "Available layer_condition_name values are:" << std::endl << rLayerConditionNames << std::endl;
}
    
} // end namespace kratos
