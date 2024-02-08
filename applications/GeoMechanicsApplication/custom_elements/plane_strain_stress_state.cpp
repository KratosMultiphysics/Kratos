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
//
#include "plane_strain_stress_state.h"

namespace Kratos
{

Matrix PlaneStrainStressState::CalculateBMatrix(const Matrix&         GradNpT,
                                                const Vector&         Np,
                                                const Geometry<Node>& rGeometry) const
{
    const auto dimension       = rGeometry.WorkingSpaceDimension();
    const auto number_of_nodes = rGeometry.size();
    Matrix     result = ZeroMatrix(VOIGT_SIZE_2D_AXISYMMETRIC, dimension * number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        auto index = dimension * i;

        result(INDEX_2D_PLANE_STRAIN_XX, index + INDEX_X) = GradNpT(i, INDEX_X);
        result(INDEX_2D_PLANE_STRAIN_YY, index + INDEX_Y) = GradNpT(i, INDEX_Y);
        result(INDEX_2D_PLANE_STRAIN_XY, index + INDEX_X) = GradNpT(i, INDEX_Y);
        result(INDEX_2D_PLANE_STRAIN_XY, index + INDEX_Y) = GradNpT(i, INDEX_X);
    }

    return result;
}

double PlaneStrainStressState::CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                               double detJ,
                                                               const Geometry<Node>&) const
{
    return rIntegrationPoint.Weight() * detJ;
}

unique_ptr<StressStatePolicy> PlaneStrainStressState::Clone() const
{
    return std::make_unique<PlaneStrainStressState>();
}

Vector PlaneStrainStressState::CalculateGreenLagrangeStrain(const Matrix& rTotalDeformationGradient) const
{
    return Kratos::Vector();
}

} // namespace Kratos
