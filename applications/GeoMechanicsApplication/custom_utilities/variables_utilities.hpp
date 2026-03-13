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
    template <typename NodeContainerType, typename DataType, typename OutputIt>
    static OutputIt GetNodalValues(const NodeContainerType& rNodes, const Variable<DataType>& rNodalVariable, OutputIt FirstOut)
    {
        return std::ranges::transform(rNodes, FirstOut, [&rNodalVariable](const auto& rNode) {
            return rNode.FastGetSolutionStepValue(rNodalVariable);
        }).out;
    }

    template <unsigned int TNumNodes, typename NodeContainerType>
    static array_1d<double, TNumNodes> GetNodalValues(const NodeContainerType& rNodes,
                                                      const Variable<double>&  rNodalVariable)
    {
        auto result = array_1d<double, TNumNodes>{};
        GetNodalValues(rNodes, rNodalVariable, result.begin());
        return result;
    }

    template <typename NodeContainerType, typename DataType>
    static std::vector<DataType> GetNodalValues(const NodeContainerType& rNodes, const Variable<DataType>& rVariable)
    {
        auto result = std::vector<DataType>{};
        result.reserve(rNodes.size());
        GetNodalValues(rNodes, rVariable, std::back_inserter(result));
        return result;
    }

    static const Variable<double>& GetComponentFromVectorVariable(const std::string& rSourceVariableName,
                                                                  const std::string& rComponent);
};

} // namespace Kratos
