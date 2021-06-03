//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes
      
// External includes
      
// Project includes
#include "mdpa_io_modeler.h"
#include "includes/model_part_io.h"


namespace Kratos
{
    MdpaIoModeler::MdpaIoModeler(Model& rModel, Parameters ModelerParameters)
        : Modeler(rModel, ModelerParameters)
        , mpModel(&rModel)
    {
        mParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    }

    const Parameters MdpaIoModeler::GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "echo_level"                                 : 0,
            "input_file_name"                            : "",
            "model_part_name"                            : "model_part",
            "skip_timer"                                 : true,
            "ignore_variables_not_in_solution_step_data" : false
        })");
        return default_parameters;
    }

    void MdpaIoModeler::SetupGeometryModel()
    {
        const auto input_file_name = mParameters["input_file_name"].GetString();
        const auto model_part_name = mParameters["model_part_name"].GetString();
        ModelPart& model_part =
            mpModel->HasModelPart(model_part_name) ?
            mpModel->GetModelPart(model_part_name) :
            mpModel->CreateModelPart(model_part_name);
        Flags options = IO::READ;
        if (mParameters["skip_timer"].GetBool())
            options = IO::SKIP_TIMER | options;
        if (mParameters["ignore_variables_not_in_solution_step_data"].GetBool())
            options = IO::IGNORE_VARIABLES_ERROR | options;
        KRATOS_INFO_IF("::[MdpaIoModeler]::", mParameters["echo_level"].GetInt() > 0) << "Importing mdpa from: " << input_file_name << std::endl;
        ModelPartIO(input_file_name, options).ReadModelPart(model_part);
    }
}
