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

// Project includes
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

#include <string>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) ProcessUtilities
{
public:
    static std::vector<std::reference_wrapper<ModelPart>> GetModelPartsFromSettings(
        Model& rModel, const Parameters& rProcessSettings, const std::string& rProcessInfo);
};
}; // namespace Kratos
