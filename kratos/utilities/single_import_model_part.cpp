//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <filesystem>

// External includes

// Project includes
#include "includes/model_part_io.h"
#include "includes/reorder_consecutive_model_part_io.h"
#include "utilities/single_import_model_part.h"
#include "processes/reorder_and_optimize_modelpart_process.h"

namespace Kratos
{
void SingleImportModelPart::Import(
    ModelPart& rModelPart,
    Parameters ModelPartImportParameters,
    const std::string& InputType
    )
{
    // Check type file
    if (InputType == "mdpa") {
        auto default_settings = Parameters(R"({
            "input_filename"                             : "",
            "skip_timer"                                 : true,
            "ignore_variables_not_in_solution_step_data" : false,
            "reorder"                                    : false,
            "reorder_consecutive"                        : false
        })");

        // Cannot validate as this might contain other settings too
        ModelPartImportParameters.AddMissingParameters(default_settings);

        const auto& input_filename = ModelPartImportParameters["input_filename"].GetString();

        // Setting some mdpa-import-related flags
        auto import_flags = ModelPartIO::READ;

        if (ModelPartImportParameters["skip_timer"].GetBool()) {
            import_flags = ModelPartIO::SKIP_TIMER|import_flags;
        }

        if (ModelPartImportParameters["ignore_variables_not_in_solution_step_data"].GetBool()) {
            import_flags = ModelPartIO::IGNORE_VARIABLES_ERROR|import_flags;
        }

        KRATOS_INFO("SingleImportModelPart") << "Reading model part from file: " << std::filesystem::current_path() / input_filename << ".mdpa" << std::endl;

        if (ModelPartImportParameters["reorder_consecutive"].GetBool()) {
            ReorderConsecutiveModelPartIO(input_filename, import_flags).ReadModelPart(rModelPart);
        } else {
            ModelPartIO(input_filename, import_flags).ReadModelPart(rModelPart);
        }

        if (ModelPartImportParameters["reorder"].GetBool()) {
            auto tmp = Parameters(R"({})");
            ReorderAndOptimizeModelPartProcess(rModelPart, tmp).Execute();
        }

        KRATOS_INFO("SingleImportModelPart") << "Finished reading model part from mdpa file." << std::endl;
    } else {
        KRATOS_ERROR << "Other model part input options are not yet implemented. Demanded: " << InputType << std::endl;
    }
}

}