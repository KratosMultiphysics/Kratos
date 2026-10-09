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
//                   Markelov Gennady
//

#pragma once

// Project includes
#include "includes/model_part.h"
#include <string>

using namespace std::string_literals;

namespace Kratos
{

class Model;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ProcessUtilities
{
public:
    static std::vector<std::reference_wrapper<ModelPart>> GetModelPartsFromSettings(
        Model&             rModel,
        const Parameters&  rProcessSettings,
        const std::string& rProcessInfo,
        const std::vector<std::string>& rModelPartNameKeys = {"model_part_name"s, "model_part_name_list"s});

    static void AddProcessesSubModelPartListToSolverSettings(const Parameters& rProjectParameters,
                                                             Parameters&       rSolverSettings);
};
}; // namespace Kratos
