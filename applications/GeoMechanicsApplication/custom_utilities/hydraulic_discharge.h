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

#include "custom_utilities/node_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/kratos_export_api.h"
#include "includes/properties.h"

namespace Kratos
{
class KRATOS_API(GEO_MECHANICS_APPLICATION) HydraulicDischarge
{
public:
    static void CalculateHydraulicDischarge(const std::vector<array_1d<double, 3>>& rFluidFlux,
                                            const std::vector<double>& rIntegrationCoefficients,
                                            const Geometry<Node>::ShapeFunctionsGradientsType& rdNDxContainer,
                                            GeometryData::IntegrationMethod IntegrationMethod,
                                            Geometry<Node>&                 rGeometry);
};

} // namespace Kratos