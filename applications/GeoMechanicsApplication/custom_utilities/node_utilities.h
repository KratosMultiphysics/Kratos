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

#include "containers/array_1d.h"
#include "containers/variable.h"
#include "includes/model_part.h"

#include <map>

namespace Kratos
{

class Node;

class KRATOS_API(GEO_MECHANICS_APPLICATION) NodeUtilities
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(NodeUtilities);

    static void AssignUpdatedVectorVariableToNonFixedComponents(Node& rNode,
                                                                const Variable<array_1d<double, 3>>& rDestinationVariable,
                                                                const array_1d<double, 3>& rNewValues,
                                                                IndexType SolutionStepIndex = 0);

    static void AssignUpdatedVectorVariableToNodes(const ModelPart::NodesContainerType& rNodes,
                                                   const Variable<array_1d<double, 3>>& rDestinationVariable,
                                                   const array_1d<double, 3>& rNewValues,
                                                   IndexType SolutionStepIndex = 0);

    static std::map<IndexType, IndexType> CreateGlobalToLocalNodeIndexMap(const PointerVector<Node>& rNodes);

    template <class TValueType>
    static void ThreadSafeNodeWrite(Node& rNode, const Variable<TValueType>& Var, const TValueType Value)
    {
        rNode.SetLock();
        rNode.FastGetSolutionStepValue(Var) = Value;
        rNode.UnSetLock();
    }
};

} // namespace Kratos
