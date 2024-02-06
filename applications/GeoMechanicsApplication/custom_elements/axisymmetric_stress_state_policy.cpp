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

#include "axisymmetric_stress_state_policy.h"
#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{

static const SizeType NumberOfUDofPerNode = 2; // Axi-symmetric stress state assumes a 2D problem definition

Matrix AxisymmetricStressStatePolicy::CalculateBMatrix(const Matrix&         GradNpT,
                                                       const Vector&         Np,
                                                       const Geometry<Node>& rGeometry) const
{
    const double radius = GeoElementUtilities::CalculateRadius(Np, rGeometry);

    const auto dimension       = rGeometry.WorkingSpaceDimension();
    const auto number_of_nodes = rGeometry.size();
    Matrix result = ZeroMatrix(VOIGT_SIZE_2D_AXISYMMETRIC, NumberOfUDofPerNode * number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const IndexType index = dimension * i;

        result(INDEX_2D_PLANE_STRAIN_XX, index + INDEX_X) = GradNpT(i, INDEX_X);
        result(INDEX_2D_PLANE_STRAIN_YY, index + INDEX_Y) = GradNpT(i, INDEX_Y);
        result(INDEX_2D_PLANE_STRAIN_ZZ, index + INDEX_X) = Np[i] / radius;
        result(INDEX_2D_PLANE_STRAIN_XY, index + INDEX_X) = GradNpT(i, INDEX_Y);
        result(INDEX_2D_PLANE_STRAIN_XY, index + INDEX_Y) = GradNpT(i, INDEX_X);
    }

    return result;
}

double AxisymmetricStressStatePolicy::CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                                      double detJ) const
{
    return 0;
}

} // namespace Kratos
