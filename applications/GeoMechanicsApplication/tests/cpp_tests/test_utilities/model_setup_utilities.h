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

namespace Kratos
{
class Model;
class ModelPart;
} // namespace Kratos

namespace Kratos::Testing::ModelSetupUtilities
{

ModelPart& CreateModelPartWithASingle2D3NElement(Model& rModel);
ModelPart& CreateModelPartWithASingle3D4NElement(Model& rModel);

} // namespace Kratos::ModelSetupUtilities
