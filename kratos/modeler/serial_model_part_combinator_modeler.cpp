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

// Project includes
#include "includes/define.h"
#include "modeler/serial_model_part_combinator_modeler.h"
#include "utilities/single_import_model_part.h"
#include "utilities/model_part_combination_utilities.h"

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

    // Multiple import
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
        SingleImportModelPart::Import(r_aux_model_part, r_model_import_settings, input_type);
    }
    ModelPartCombinationUtilities(*mpModel).CombineModelParts(combine_param);
}

} // namespace Kratos
