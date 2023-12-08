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

#include "includes/ublas_interface.h"

namespace Kratos
{

class VariablesUtilities
{
public:
    template <typename GeometryType>
    static Vector GetNodalValues(const GeometryType& rGeometry,
                                 const Variable<double>& rNodalVariable)
    {
        Vector result{rGeometry.size(), 0.0};
        std::transform(rGeometry.begin(), rGeometry.end(), result.begin(),
                       [&rNodalVariable](const auto& node) {
                           return node.FastGetSolutionStepValue(rNodalVariable);
                       });
        return result;
    }
};

}
