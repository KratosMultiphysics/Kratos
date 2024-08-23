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

#include "integration/integration_point.h"

#include <vector>

namespace Kratos
{

class IntegrationScheme
{
public:
    using IntegrationPointType       = IntegrationPoint<3>;
    using IntegrationPointVectorType = std::vector<IntegrationPointType>;
};

} // namespace Kratos
