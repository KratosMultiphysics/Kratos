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

InputUtilityStub::InputUtilityStub(const std::string& parameter_json_string) :
    mParameterJsonString{parameter_json_string}
{

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

int InputUtilityStub::NumberOfReadCalls() const {
    return mNumberOfReadCalls;
}

int InputUtilityStub::NumberOfMaterialCalls() const {
    return mNumberOfMaterialCalls;
}

} // Kratos