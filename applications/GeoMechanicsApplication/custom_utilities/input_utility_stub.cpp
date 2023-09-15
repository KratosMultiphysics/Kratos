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

InputUtilityStub::InputUtilityStub()
        : mNumberOfReadCalls{0},
        mNumberOfMaterialCalls{0}
{
    mParameterJsonString = "{"
                           "\"solver_settings\":"
                           "{"
                           "\"model_part_name\":\"test\","
                           "\"model_import_settings\":"
                           "{"
                           "\"input_type\": \"mdpa\","
                           "\"input_filename\": \"mesh_stage1\""
                           "},"
                           "\"material_import_settings\": "
                           "{"
                           "\"materials_filename\": \"MaterialParameters1.json\""
                           "}"
                           "}"
                           "}";
}

Parameters InputUtilityStub::ProjectParametersFrom(const std::string& rProjectFilePath) const
{
    Parameters result(mParameterJsonString);

    return result;
}

void InputUtilityStub::ReadModelFromFile(const std::filesystem::path& rModelPartFilePath, ModelPart& rModelPart) const
{
    mNumberOfReadCalls++;
}

void InputUtilityStub::AddMaterialsFrom(const std::string& rMaterialFilePath, Model& rModel) const
{
    mNumberOfMaterialCalls++;
}

int InputUtilityStub::numberOfReadCalls() const {
    return mNumberOfReadCalls;
}

int InputUtilityStub::numberOfMaterialCalls() const {
    return mNumberOfMaterialCalls;
}

} // Kratos