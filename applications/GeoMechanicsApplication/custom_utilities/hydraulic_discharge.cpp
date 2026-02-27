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

#include "hydraulic_discharge.h"
#include "custom_utilities/node_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
void HydraulicDischarge::CalculateHydraulicDischarge(const std::vector<array_1d<double, 3>>& rFluidFlux,
                                                     const std::vector<double>& rIntegrationCoefficients,
                                                     const Geometry<Node>::ShapeFunctionsGradientsType& rdNDxContainer,
                                                     GeometryData::IntegrationMethod IntegrationMethod,
                                                     Geometry<Node>& rGeometry)
{
    KRATOS_TRY

    const IndexType number_of_integration_points = rGeometry.IntegrationPointsNumber(IntegrationMethod);
    const auto dimension        = rGeometry.WorkingSpaceDimension();
    const auto number_of_points = rGeometry.PointsNumber();
    Matrix     grad_Np_T(number_of_points, dimension);

    for (unsigned int integration_point = 0; integration_point < number_of_integration_points; ++integration_point) {
        noalias(grad_Np_T) = rdNDxContainer[integration_point];

        const auto integration_coefficient = rIntegrationCoefficients[integration_point];

        for (unsigned int node = 0; node < number_of_points; ++node) {
            double hydraulic_discharge = 0;
            for (unsigned int direction = 0; direction < dimension; ++direction) {
                hydraulic_discharge += grad_Np_T(node, direction) * rFluidFlux[integration_point][direction];
            }

            hydraulic_discharge *= integration_coefficient;
            hydraulic_discharge += rGeometry[node].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE);
            NodeUtilities::ThreadSafeNodeWrite(rGeometry[node], HYDRAULIC_DISCHARGE, hydraulic_discharge);
        }
    }

    KRATOS_CATCH("")
}

} // namespace Kratos