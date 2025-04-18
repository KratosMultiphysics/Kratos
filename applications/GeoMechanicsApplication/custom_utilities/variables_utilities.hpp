// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Anne van de Graaf
//

#pragma once

#include <algorithm>

#include "containers/variable.h"
#include "geometries/geometry.h"
#include "includes/node.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) VariablesUtilities
{
public:
    template <typename OutputIt>
    static OutputIt GetNodalValues(const Geometry<Node>& rGeometry, const Variable<double>& rNodalVariable, OutputIt FirstOut)
    {
        return std::transform(rGeometry.begin(), rGeometry.end(), FirstOut, [&rNodalVariable](const auto& node) {
            return node.FastGetSolutionStepValue(rNodalVariable);
        });
    }

    static const Variable<double>& GetComponentFromVectorVariable(const std::string& rSourceVariableName,
                                                                  const std::string& rComponent);
};

} // namespace Kratos
