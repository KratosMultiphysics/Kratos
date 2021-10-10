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

// Project includes
#include "includes/define.h"
#include "utilities/model_part_combination_utilities.h"
#include "modeler/serial_model_part_combinator_modeler.h"

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
    // if input_type.IsArray():
    //     if model_part_import_settings["input_filename"].IsArray():
    //         current_model = model_part.GetModel()
    //         combine_param = KratosMultiphysics.Parameters("""{
    //             "model_parts_list"         : [],
    //             "combined_model_part_name" : ""
    //         }""")
    //         combine_param["combined_model_part_name"].SetString(model_part.Name)
    //         filenames_list = model_part_import_settings["input_filename"].GetStringArray()
    //         copy_model_part_import_settings = KratosMultiphysics.Parameters(model_part_import_settings)
    //         copy_model_part_import_settings.RemoveValue("input_filename")
    //         copy_model_part_import_settings.AddString("input_filename", "")
    //         for i in range(len(filenames_list)):
    //             aux_name = "AUX_MODELPART" + str(i)
    //             combine_param["model_parts_list"].Append(aux_name)
    //             aux_model_part = current_model.CreateModelPart(aux_name)
    //             copy_model_part_import_settings["input_filename"].SetString(filenames_list[i])
    //             self._single_ImportModelPart(aux_model_part, model_part_import_settings, input_type[i].GetString())
    //         KratosMultiphysics.ModelPartCombinationUtilities(current_model).CombineModelParts(combine_param)
    //     else:
    //         raise Exception("Multiple input formats, but only one file?.")
    // elif input_type.GetString() == "mdpa": # Add other file reading formats
    //     if model_part_import_settings["input_filename"].IsArray():
    //         current_model = model_part.GetModel()
    //         combine_param = KratosMultiphysics.Parameters("""{
    //             "model_parts_list"         : [],
    //             "combined_model_part_name" : ""
    //         }""")
    //         combine_param["combined_model_part_name"].SetString(model_part.Name)
    //         filenames_list = model_part_import_settings["input_filename"].GetStringArray()
    //         copy_model_part_import_settings = KratosMultiphysics.Parameters(model_part_import_settings)
    //         copy_model_part_import_settings.RemoveValue("input_filename")
    //         copy_model_part_import_settings.AddString("input_filename", "")
    //         for i in range(len(filenames_list)):
    //             aux_name = "AUX_MODELPART" + str(i)
    //             combine_param["model_parts_list"].Append(aux_name)
    //             aux_model_part = current_model.CreateModelPart(aux_name)
    //             copy_model_part_import_settings["input_filename"].SetString(filenames_list[i])
    //             self._single_ImportModelPart(aux_model_part, model_part_import_settings, input_type.GetString())
    //         KratosMultiphysics.ModelPartCombinationUtilities(current_model).CombineModelParts(combine_param)
    //     else:
    //         self._single_ImportModelPart(model_part, model_part_import_settings, input_type.GetString())

}

} // namespace Kratos
