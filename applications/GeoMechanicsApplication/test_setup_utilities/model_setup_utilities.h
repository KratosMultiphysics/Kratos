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

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "includes/kratos_export_api.h"

namespace Kratos
{
class Model;
class ModelPart;
} // namespace Kratos

namespace Kratos::Testing
{
class KRATOS_API(KRATOS_GEO_TEST_UTILS) ModelSetupUtilities
{
public:
    static ModelPart& CreateModelPartWithASingle2D3NElement(Model& rModel,
                                                            const Geo::ConstVariableRefs& rNodalVariables = {});

    static ModelPart& CreateModelPartWithASingle2D6NDiffOrderElement(Model& rModel);

    static ModelPart& CreateModelPartWithASingle3D4NElement(Model& rModel,
                                                            const Geo::ConstVariableRefs& rNodalVariables = {});

    static ModelPart& CreateModelPartWithASingle3D8NElement(Model& rModel,
                                                            const Geo::ConstVariableRefs& rNodalVariables = {});

    static ModelPart& CreateModelPartWithASingle3D20NElement(Model& rModel,
                                                             const Geo::ConstVariableRefs& rNodalVariables = {});

    static ModelPart& CreateModelPartWithASingle2D2NElement(Model& rModel,
                                                            const Geo::ConstVariableRefs& rNodalVariables = {});

    static ModelPart& CreateModelPartWithASingle2D10NElement(Model& rModel,
                                                             const Geo::ConstVariableRefs& rNodalVariables = {});

    static ModelPart& CreateModelPartWithASingle2D15NElement(Model& rModel,
                                                             const Geo::ConstVariableRefs& rNodalVariables = {});

    static ModelPart& CreateModelPartWithASingle3D6NInterfaceElement(
        Model& rModel, const Geo::ConstVariableRefs& rNodalVariables = {});

    static ModelPart&          CreateModelPartWithASingle2D6NUPwDiffOrderElement(Model& rModel);
    static ModelPart&          CreateModelPartWithASingle3D10NUPwDiffOrderElement(Model& rModel);
    static Triangle2D3<Node>   Create2D3NTriangleGeometry();
    static Tetrahedra3D4<Node> Create3D4NTetrahedraGeometry();
    static void CreateNumberOfNewNodes(ModelPart& rModelPart, std::size_t NumberOfNodes);
    static PointerVector<Node> GetNodesFromIds(ModelPart& rModelPart, const std::vector<std::size_t>& rNodeIds);
};

} // namespace Kratos::Testing
