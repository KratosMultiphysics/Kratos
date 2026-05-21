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

#pragma once

#include "custom_utilities/input_utility.h"

namespace Kratos
{

class StubInputUtility : public InputUtility
{
public:
    explicit StubInputUtility(const std::string& rParameterJsonString);

    [[nodiscard]] Parameters ProjectParametersFromFile(const std::filesystem::path& rProjectFilePath) const override;
    void ReadModelFromFile(const std::filesystem::path& rModelPartFilePath, ModelPart& rModelPart) const override;
    void AddMaterialsFromFile(const std::filesystem::path& rMaterialFilePath, Model& rModel) const override;

    int NumberOfReadCalls() const;
    int NumberOfMaterialCalls() const;

private:
    mutable int mNumberOfReadCalls     = 0;
    mutable int mNumberOfMaterialCalls = 0;
    std::string mParameterJsonString;
};

} // namespace Kratos