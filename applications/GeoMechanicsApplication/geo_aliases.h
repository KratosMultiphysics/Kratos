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
#include <memory>
#include <optional>
#include <vector>

#include "containers/variable.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "integration/integration_point.h"

namespace Kratos::Geo
{

using ConstVariableRefs     = std::vector<std::reference_wrapper<const Variable<double>>>;
using ConstVariableDataRefs = std::vector<std::reference_wrapper<const VariableData>>;

using IntegrationPointType       = IntegrationPoint<3>;
using IntegrationPointVectorType = std::vector<IntegrationPointType>;

using KappaDependentFunction = std::function<double(double)>;

using GeometryUniquePtr         = std::unique_ptr<Geometry<Node>>;
using OptionalGeometryUniquePtr = std::optional<GeometryUniquePtr>;

} // namespace Kratos::Geo
