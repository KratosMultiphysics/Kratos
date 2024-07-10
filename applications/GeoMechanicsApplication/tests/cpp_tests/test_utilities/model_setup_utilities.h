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

#include <geometries/tetrahedra_3d_4.h>
#include <geometries/triangle_2d_3.h>

namespace Kratos
{
class Model;
class ModelPart;
} // namespace Kratos

namespace Kratos::Testing::ModelSetupUtilities
{

ModelPart& CreateModelPartWithASingle2D3NElement(Model& rModel,
                                                 const Geo::ConstVariableRefs& rNodalVariables = {});

ModelPart& CreateModelPartWithASingle2D6NDiffOrderElement(Model& rModel);

ModelPart& CreateModelPartWithASingle3D4NElement(Model& rModel,
                                                 const Geo::ConstVariableRefs& rNodalVariables = {});

ModelPart& CreateModelPartWithASingle2D6NUPwDiffOrderElement(Model& rModel);
ModelPart& CreateModelPartWithASingle3D10NUPwDiffOrderElement(Model& rModel);
Triangle2D3<Node> Create2D3NTriangleGeometry();
Tetrahedra3D4<Node> Create3D4NTetrahedraGeometry();

} // namespace Kratos::Testing::ModelSetupUtilities
