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

    template <typename NodeContainerType, typename DataType>
    static std::vector<DataType> GetNodalValues(const NodeContainerType& rNodes, const Variable<DataType>& rVariable)
    {
        auto result = std::vector<DataType>{};
        result.reserve(rNodes.size());
        auto get_nodal_value = [&rVariable](const auto& rNode) {
            return rNode.FastGetSolutionStepValue(rVariable);
        };
        std::ranges::transform(rNodes, std::back_inserter(result), get_nodal_value);
        return result;
    }

    static const Variable<double>& GetComponentFromVectorVariable(const std::string& rSourceVariableName,
                                                                  const std::string& rComponent);
};

} // namespace Kratos
