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

#include <filesystem>
#include <includes/kratos_parameters.h>
#include "containers/model.h"

namespace Kratos {

class InputUtility {
public:
    virtual ~InputUtility() = default;

    [[nodiscard]] virtual Parameters ProjectParametersFromFile(const std::filesystem::path& rProjectFilePath) const = 0;
    virtual void ReadModelFromFile(const std::filesystem::path& rModelPartFilePath, ModelPart& rModelPart) const = 0;
    virtual void AddMaterialsFromFile(const std::filesystem::path& rMaterialFilePath, Model& rModel) const = 0;
};

}
