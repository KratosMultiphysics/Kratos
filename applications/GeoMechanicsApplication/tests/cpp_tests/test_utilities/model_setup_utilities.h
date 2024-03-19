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

#include "geo_aliases.h"

namespace Kratos
{
class Model;
class ModelPart;
} // namespace Kratos

namespace Kratos::Testing::ModelSetupUtilities
{

ModelPart& CreateModelPartWithASingle2D3NElement(Model& rModel,
                                                 const Geo::ConstVariableRefs& rNodalVariables = {});
ModelPart& CreateModelPartWithASingle3D4NElement(Model& rModel,
                                                 const Geo::ConstVariableRefs& rNodalVariables = {});

ModelPart& CreateModelPartWithASingle2D6NUPwDiffOrderElement(Model& rModel);
ModelPart& CreateModelPartWithASingle3D10NUPwDiffOrderElement(Model& rModel);

} // namespace Kratos::Testing::ModelSetupUtilities
