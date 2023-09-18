// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//                   Richard Faasse
//

#pragma once

#include "interface_input_utility.h"

namespace Kratos {

class InputUtility : public InterfaceInputUtility {
public:
    [[nodiscard]] Parameters ProjectParametersFrom(const std::string& rProjectFilePath) const override;
    void ReadModelFromFile(const std::filesystem::path& rModelPartFilePath, ModelPart& rModelPart) const override;
    void AddMaterialsFrom(const std::string& rMaterialFilePath, Model& rModel) const override;
};

}
