// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#pragma once

#include <functional>
#include <vector>

#include "containers/variable.h"
#include "integration/integration_point.h"

namespace Kratos::Geo
{

using ConstVariableRefs = std::vector<std::reference_wrapper<const Variable<double>>>;

using IntegrationPointType       = IntegrationPoint<3>;
using IntegrationPointVectorType = std::vector<IntegrationPointType>;

} // namespace Kratos::Geo
