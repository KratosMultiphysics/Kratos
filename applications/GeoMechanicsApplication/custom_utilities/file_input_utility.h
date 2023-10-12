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

#include "input_utility.h"

namespace Kratos {

class FileInputUtility : public InputUtility {
public:
    [[nodiscard]] Parameters ProjectParametersFromFile(const std::filesystem::path& rProjectFilePath) const override;
    void ReadModelFromFile(const std::filesystem::path& rModelPartFilePath, ModelPart& rModelPart) const override;
    void AddMaterialsFromFile(const std::filesystem::path& rMaterialFilePath, Model& rModel) const override;
};

}
