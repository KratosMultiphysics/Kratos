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
#include <vector>

#include "containers/variable.h"
#include "includes/constitutive_law.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "integration/integration_point.h"

namespace Kratos
{
class ProcessInfo;
class Properties;
} // namespace Kratos

namespace Kratos::Geo
{

using ConstVariableRefs     = std::vector<std::reference_wrapper<const Variable<double>>>;
using ConstVariableDataRefs = std::vector<std::reference_wrapper<const VariableData>>;

using IntegrationPointType       = IntegrationPoint<3>;
using IntegrationPointVectorType = std::vector<IntegrationPointType>;

using KappaDependentFunction = std::function<double(double)>;

using GeometryUniquePtr = std::unique_ptr<Geometry<Node>>;

using BMatrixVectorGetter           = std::function<std::vector<Matrix>()>;
using StrainVectorsGetter           = std::function<std::vector<Vector>()>;
using IntegrationCoefficientsGetter = std::function<std::vector<double>()>;
using PropertiesGetter              = std::function<const Properties&()>;
using ProcessInfoGetter             = std::function<const ProcessInfo&()>;
using ConstitutiveLawVectorGetter   = std::function<const std::vector<ConstitutiveLaw::Pointer>&()>;
} // namespace Kratos::Geo
