// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Marjan Fathian
//

#include "axisymmetric_stress_state.h"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

Matrix AxisymmetricStressState::CalculateBMatrix(const Matrix&         rGradNpT,
                                                 const Vector&         rNp,
                                                 const Geometry<Node>& rGeometry) const
{
    const auto radius = GeoElementUtilities::CalculateRadius(rNp, rGeometry);

    const auto dimension       = rGeometry.WorkingSpaceDimension();
    const auto number_of_nodes = rGeometry.size();
    Matrix     result = ZeroMatrix(VOIGT_SIZE_2D_AXISYMMETRIC, dimension * number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const IndexType index = dimension * i;

        result(INDEX_2D_PLANE_STRAIN_XX, index + INDEX_X) = rGradNpT(i, INDEX_X);
        result(INDEX_2D_PLANE_STRAIN_YY, index + INDEX_Y) = rGradNpT(i, INDEX_Y);
        result(INDEX_2D_PLANE_STRAIN_ZZ, index + INDEX_X) = rNp[i] / radius;
        result(INDEX_2D_PLANE_STRAIN_XY, index + INDEX_X) = rGradNpT(i, INDEX_Y);
        result(INDEX_2D_PLANE_STRAIN_XY, index + INDEX_Y) = rGradNpT(i, INDEX_X);
    }

    return result;
}

double AxisymmetricStressState::CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                                double DetJ,
                                                                const Geometry<Node>& rGeometry) const
{
    Vector shape_function_values;
    shape_function_values =
        rGeometry.ShapeFunctionsValues(shape_function_values, rIntegrationPoint.Coordinates());

    const auto radius_weight =
        GeoElementUtilities::CalculateAxisymmetricCircumference(shape_function_values, rGeometry);

    return rIntegrationPoint.Weight() * DetJ * radius_weight;
}

std::unique_ptr<StressStatePolicy> AxisymmetricStressState::Clone() const
{
    return std::make_unique<AxisymmetricStressState>();
}

Vector AxisymmetricStressState::CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const
{
    KRATOS_ERROR << "The calculation of Green Lagrange strain is not implemented for axisymmetric "
                    "configurations.\n";
}

} // namespace Kratos
