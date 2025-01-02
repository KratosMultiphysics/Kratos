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
#include "cad_io_thb_modeler.h"
// #include "input_output/cad_json_input.h"
#include "input_output/cad_json_thb_input.h"
#include "input_output/cad_json_output.h"


namespace Kratos
{
    ///@name Stages
    ///@{

    void CadIoThbModeler::SetupGeometryModel()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("cad_model_part_name"))
            << "Missing \"cad_model_part_name\" in CadIoThbModeler Parameters." << std::endl;
        const std::string cad_model_part_name = mParameters["cad_model_part_name"].GetString();
        ModelPart& cad_model_part = mpModel->HasModelPart(cad_model_part_name)
            ? mpModel->GetModelPart(cad_model_part_name)
            : mpModel->CreateModelPart(cad_model_part_name);

        const std::string DataFileName = mParameters.Has("geometry_file_name")
            ? mParameters["geometry_file_name"].GetString()
            : "geometry.cad.json";

        KRATOS_INFO_IF("::[CadIoThbModeler]::", mEchoLevel > 0) << "Importing Cad Model from: " << DataFileName << std::endl;

        CadJsonThbInput<Node, Point>(
        // CadJsonInput<Node, Point>(
            DataFileName, mEchoLevel).ReadModelPart(cad_model_part);
    }

    void CadIoThbModeler::SetupModelPart()
    {
        // if (mParameters.Has("output_geometry_file_name")) {
        //     std::string DataFileName = mParameters["output_geometry_file_name"].GetString();

        //     const std::string cad_model_part_name = mParameters["cad_model_part_name"].GetString();
        //     ModelPart& cad_model_part = mpModel->HasModelPart(cad_model_part_name)
        //         ? mpModel->GetModelPart(cad_model_part_name)
        //         : mpModel->CreateModelPart(cad_model_part_name);

        //     std::string output_file_text;
        //     CadJsonOutput::GetCadJsonOutput(cad_model_part, output_file_text, mEchoLevel);

        //     std::ofstream output_file(DataFileName);
        //     output_file << output_file_text;
        //     output_file.close();
        // }
    }
    ///@}
}
