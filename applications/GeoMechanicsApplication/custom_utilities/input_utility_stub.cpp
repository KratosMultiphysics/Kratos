// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//


#include "input_utility_stub.h"

namespace Kratos {

Parameters InputUtilityStub::ProjectParametersFrom(const std::string& rProjectFilePath) const
{
    const std::string jsonString = "{"
                                   "\"solver_settings\":"
                                   "{"
                                   "\"model_part_name\":\"test\","
                                   "\"model_import_settings\":"
                                   "{"
                                   "\"input_type\": \"mdpa\","
                                   "\"input_filename\": \"mesh_stage1\""
                                   "}"
                                   "}"
                                   "}";
    Parameters result(jsonString);

    return result;
}

void InputUtilityStub::ReadModelFromFile(const std::filesystem::path& rModelPartFilePath, ModelPart& rModelPart) const
{

}

void InputUtilityStub::AddMaterialsFrom(const std::string& rMaterialFilePath, Model& rModel) const
{

}

} // Kratos