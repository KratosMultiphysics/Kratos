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

namespace Kratos
{

class Model;
class Parameters;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ProcessUtilities
{
public:
    static std::vector<std::reference_wrapper<ModelPart>> GetModelPartsFromSettings(
        Model& rModel, const Parameters& rProcessSettings, const std::string& rProcessInfo);

    static void AddProcessesSubModelPartList(const Parameters& rProjectParameters, Parameters& rSolverSettings);
};
}; // namespace Kratos
