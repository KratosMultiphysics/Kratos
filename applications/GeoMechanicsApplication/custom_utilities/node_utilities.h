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

namespace Kratos
{

class Node;

class KRATOS_API(GEO_MECHANICS_APPLICATION) NodeUtilities
{
public:
    static void AssignUpdatedVectorVariableToNonFixedComponents(Node& rNode,
                                                                const Variable<array_1d<double, 3>>& rDestinationVariable,
                                                                const array_1d<double, 3>& rNewValues);
};

} // namespace Kratos
