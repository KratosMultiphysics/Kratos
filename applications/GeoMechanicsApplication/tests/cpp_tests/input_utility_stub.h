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

#include "custom_utilities/interface_input_utility.h"

namespace Kratos {

class InputUtilityStub : public InterfaceInputUtility
{
public:
    explicit InputUtilityStub(const std::string& parameter_json_string);

    [[nodiscard]] Parameters ProjectParametersFrom(const std::string& rProjectFilePath) const override;
    void ReadModelFromFile(const std::filesystem::path& rModelPartFilePath, ModelPart& rModelPart) const override;
    void AddMaterialsFrom(const std::string& rMaterialFilePath, Model& rModel) const override;

    int NumberOfReadCalls() const;
    int NumberOfMaterialCalls() const;

private:
    mutable int mNumberOfReadCalls = 0;
    mutable int mNumberOfMaterialCalls = 0;
    std::string mParameterJsonString;
};

} // Kratos
