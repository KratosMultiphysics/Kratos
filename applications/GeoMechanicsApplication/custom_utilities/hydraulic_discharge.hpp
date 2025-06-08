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

#include "includes/kratos_export_api.h"
#include "includes/properties.h"
#include "geo_mechanics_application_variables.h"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{
template<unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) HydraulicDischarge
{
public:
    static void CalculateHydraulicDischarge(const std::vector<array_1d<double, 3>>& rFluidFlux,
                                     const std::vector<double>& rIntegrationCoefficients,
                                     const Geometry<Node>::ShapeFunctionsGradientsType& rdNDxContainer,
                                     GeometryData::IntegrationMethod IntegrationMethod,
                                     Geometry<Node>& rGeometry)
    {
        KRATOS_TRY

        const IndexType     number_of_integration_points =
            rGeometry.IntegrationPointsNumber(IntegrationMethod);
        Matrix grad_Np_T(TNumNodes, TDim);

        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            noalias(grad_Np_T) = rdNDxContainer[integration_point];

            auto integration_coefficient = rIntegrationCoefficients[integration_point];

            for (unsigned int node = 0; node < TNumNodes; ++node) {
                double hydraulic_discharge = 0;
                for (unsigned int direction = 0; direction < TDim; ++direction) {
                    hydraulic_discharge +=
                        grad_Np_T(node, direction) * rFluidFlux[integration_point][direction];
                }

                hydraulic_discharge *= integration_coefficient;
                hydraulic_discharge += rGeometry[node].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE);
                GeoElementUtilities::ThreadSafeNodeWrite(rGeometry[node],
                                                         HYDRAULIC_DISCHARGE, hydraulic_discharge);
            }
             }

        KRATOS_CATCH("")
    }

};

} // namespace Kratos