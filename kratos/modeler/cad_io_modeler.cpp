//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


// Project includes
#include "cad_io_modeler.h"
#include "input_output/cad_json_input.h"
#include "input_output/cad_json_output.h"
#include "utilities/nurbs_utilities/projection_nurbs_geometry_utilities.h"


namespace Kratos
{
    ///@name Stages
    ///@{

    void CadIoModeler::SetupGeometryModel()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("cad_model_part_name"))
            << "Missing \"cad_model_part_name\" in CadIoModeler Parameters." << std::endl;

        const std::string cad_model_part_name = mParameters["cad_model_part_name"].GetString();

        ModelPart& cad_model_part = mpModel->HasModelPart(cad_model_part_name)
            ? mpModel->GetModelPart(cad_model_part_name)
            : mpModel->CreateModelPart(cad_model_part_name);

        const std::string DataFileName = mParameters.Has("geometry_file_name")
            ? mParameters["geometry_file_name"].GetString()
            : "geometry.cad.json";

        const std::string LocalRefFileName = mParameters.Has("local_ref_file_name")
            ? mParameters["local_ref_file_name"].GetString()
            : "";

        KRATOS_INFO_IF("::[CadIoModeler]::", mEchoLevel > 0)
            << "Importing Cad Model from: " << DataFileName << std::endl;

        ProjectionAlgorithm projection_algorithm =
            ProjectionAlgorithm::NewtonRaphson;

        if (mParameters.Has("projection_algorithm")) {
            const std::string algorithm_name =
                mParameters["projection_algorithm"].GetString();

            if (algorithm_name == "newton_raphson") {
                projection_algorithm = ProjectionAlgorithm::NewtonRaphson;
            } else if (algorithm_name == "levenberg_marquardt") {
                projection_algorithm = ProjectionAlgorithm::LevenbergMarquardt;
            } else {
                KRATOS_ERROR
                    << "Unknown projection algorithm \""
                    << algorithm_name
                    << "\". Available options are:"
                    << "\n  - newton_raphson"
                    << "\n  - levenberg_marquardt"
                    << std::endl;
            }
        }

        if (LocalRefFileName.empty()) {
            CadJsonInput<Node, Point>(
                DataFileName, mEchoLevel, projection_algorithm).ReadModelPart(cad_model_part);
        } else {
            CadJsonInput<Node, Point>(
                DataFileName, mEchoLevel, projection_algorithm, LocalRefFileName).ReadModelPart(cad_model_part);
        }
    }

    void CadIoModeler::SetupModelPart()
    {
        if (mParameters.Has("output_geometry_file_name")) {
            std::string DataFileName = mParameters["output_geometry_file_name"].GetString();

            const std::string cad_model_part_name = mParameters["cad_model_part_name"].GetString();
            ModelPart& cad_model_part = mpModel->HasModelPart(cad_model_part_name)
                ? mpModel->GetModelPart(cad_model_part_name)
                : mpModel->CreateModelPart(cad_model_part_name);

            std::string output_file_text;
            CadJsonOutput::GetCadJsonOutput(cad_model_part, output_file_text, mEchoLevel);

            std::ofstream output_file(DataFileName);
            output_file << output_file_text;
            output_file.close();
        }
    }
    ///@}
}
