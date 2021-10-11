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

// External includes
#if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || (defined(__cplusplus) && __cplusplus >= 201703L)) && defined(__has_include)
#if __has_include(<filesystem>) && (!defined(__MAC_OS_X_VERSION_MIN_REQUIRED) || __MAC_OS_X_VERSION_MIN_REQUIRED >= 101500)
#define GHC_USE_STD_FS
#include <filesystem>
namespace fs = std::filesystem;
#endif
#endif
#ifndef GHC_USE_STD_FS
#include <ghc/filesystem.hpp>
namespace fs = ghc::filesystem;
#endif

// Project includes
#include "includes/define.h"
#include "includes/model_part_io.h"
#include "includes/kratos_filesystem.h"
#include "includes/reorder_consecutive_model_part_io.h"
#include "utilities/model_part_combination_utilities.h"
#include "modeler/serial_model_part_combinator_modeler.h"
#include "processes/reorder_and_optimize_modelpart_process.h"

namespace Kratos
{

Modeler::Pointer SerialModelPartCombinatorModeler::Create(
    Model& rModel,
    const Parameters ModelParameters
    ) const
{
    return Kratos::make_shared<SerialModelPartCombinatorModeler>(rModel, ModelParameters);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialModelPartCombinatorModeler::SetupModelPart()
{
    // Import parameters
    const auto& r_model_import_settings = mParameters["model_import_settings"];
    const auto& r_input_type = r_model_import_settings["input_type"];
    const auto& r_input_filename = r_model_import_settings["input_filename"];

    // Multiple import
    if (r_input_type.IsArray() || r_input_filename.IsArray()) {
        auto combine_param = Parameters(R"({
            "model_parts_list"         : []
        })");
        combine_param.AddValue("combined_model_part_name", mParameters["model_part_name"]);
        const auto filenames_list = r_model_import_settings["input_filename"].GetStringArray();
        auto copy_model_part_import_settings = Parameters(r_model_import_settings);
        copy_model_part_import_settings.RemoveValue("input_filename");
        copy_model_part_import_settings.AddString("input_filename", "");
        for (std::size_t i = 0; i < filenames_list.size(); ++i) {
            const std::string aux_name = "AUX_MODELPART" + std::to_string(i);
            combine_param["model_parts_list"].Append(aux_name);
            auto& r_aux_model_part = mpModel->CreateModelPart(aux_name);
            copy_model_part_import_settings["input_filename"].SetString(filenames_list[i]);
            const std::string& input_type = r_input_type.IsArray() ? r_input_type[i].GetString() : r_input_type.GetString();
            SingleImportModelPart(r_aux_model_part, r_model_import_settings, input_type);
        }
        ModelPartCombinationUtilities(*mpModel).CombineModelParts(combine_param);
    } else { // Single import
        const auto& r_model_part_name = mParameters["model_part_name"].GetString();
        auto& r_model_part = mpModel->HasModelPart(r_model_part_name) ? mpModel->GetModelPart(r_model_part_name) : mpModel->CreateModelPart(r_model_part_name);
        SingleImportModelPart(r_model_part, r_model_import_settings, r_input_type.GetString());
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SerialModelPartCombinatorModeler::SingleImportModelPart(
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

        KRATOS_INFO("SerialModelPartCombinatorModeler") << "Reading model part from file: " << FilesystemExtensions::CurrentWorkingDirectory()<< "/" << input_filename << ".mdpa" << std::endl;

        if (ModelPartImportParameters["reorder_consecutive"].GetBool()) {
            ReorderConsecutiveModelPartIO(input_filename, import_flags).ReadModelPart(rModelPart);
        } else {
            ModelPartIO(input_filename, import_flags).ReadModelPart(rModelPart);
        }

        if (ModelPartImportParameters["reorder"].GetBool()) {
            auto tmp = Parameters(R"({})");
            ReorderAndOptimizeModelPartProcess(rModelPart, tmp).Execute();
        }

        KRATOS_INFO("SerialModelPartCombinatorModeler") << "Finished reading model part from mdpa file." << std::endl;
    } else {
        KRATOS_ERROR << "Other model part input options are not yet implemented. Demanded: " << InputType << std::endl;
    }
}

} // namespace Kratos
