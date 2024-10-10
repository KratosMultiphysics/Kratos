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

#include "stub_input_utility.h"

namespace Kratos
{

StubInputUtility::StubInputUtility(const std::string& rParameterJsonString)
    : mParameterJsonString{rParameterJsonString}
{
}

Parameters StubInputUtility::ProjectParametersFromFile(const std::filesystem::path& rProjectFilePath) const
{
    return Parameters{mParameterJsonString};
}

void StubInputUtility::ReadModelFromFile(const std::filesystem::path& rModelPartFilePath, ModelPart& rModelPart) const
{
    ++mNumberOfReadCalls;
}

void StubInputUtility::AddMaterialsFromFile(const std::filesystem::path& rMaterialFilePath, Model& rModel) const
{
    ++mNumberOfMaterialCalls;
}

int StubInputUtility::NumberOfReadCalls() const { return mNumberOfReadCalls; }

int StubInputUtility::NumberOfMaterialCalls() const { return mNumberOfMaterialCalls; }

} // namespace Kratos