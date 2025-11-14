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

    template <unsigned int TNumNodes>
    static array_1d<double, TNumNodes> GetNodalValuesOf(const Variable<double>& rNodalVariable,
                                                        const Geometry<Node>&   rGeometry)
    {
        auto result = array_1d<double, TNumNodes>{};
        GetNodalValues(rGeometry, rNodalVariable, result.begin());
        return result;
    }

    static const Variable<double>& GetComponentFromVectorVariable(const std::string& rSourceVariableName,
                                                                  const std::string& rComponent);
};

} // namespace Kratos
