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

        KRATOS_INFO_IF("::[CadIoModeler]::", mEchoLevel > 0) << "Importing Cad Model from: " << DataFileName << std::endl;

        CadJsonInput<Node<3>, Point>(
            DataFileName, mEchoLevel).ReadModelPart(cad_model_part);
    }

    ///@}
}
