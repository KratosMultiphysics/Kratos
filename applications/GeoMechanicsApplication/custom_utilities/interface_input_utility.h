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

#include <string>
#include <includes/kratos_parameters.h>
#include "containers/model.h"

namespace Kratos {

class InterfaceInputUtility {
public:
    [[nodiscard]] virtual Parameters ProjectParametersFrom(const std::string &rProjectFilePath) const = 0;
    virtual void ReadModelFromFile(const std::filesystem::path& rModelPartFilePath, ModelPart& rModelPart) const = 0;
    virtual void AddMaterialsFrom(const std::string& rMaterialFilePath, Model& rModel) const = 0;
};

}
