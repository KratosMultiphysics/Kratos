// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Marjan Fathian
//                   Richard Faasse
//
#include "three_dimension_stress_state.h"

namespace Kratos
{

Matrix ThreeDimensionStressState::CalculateBMatrix(const Matrix&         GradNpT,
                                       const Vector&         Np,
                                       const Geometry<Node>& rGeometry) const
{
    const auto dimension       = rGeometry.WorkingSpaceDimension();
    const auto number_of_nodes = rGeometry.size();
    Matrix result = ZeroMatrix(VOIGT_SIZE_3D, dimension * number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
    auto index = dimension * i;

    result(INDEX_3D_XX, index + INDEX_X) = GradNpT(i, INDEX_X);
    result(INDEX_3D_YY, index + INDEX_Y) = GradNpT(i, INDEX_Y);
    result(INDEX_3D_ZZ, index + INDEX_Z) = GradNpT(i, INDEX_Z);
    result(INDEX_3D_XY, index + INDEX_X) = GradNpT(i, INDEX_Y);
    result(INDEX_3D_XY, index + INDEX_Y) = GradNpT(i, INDEX_X);
    result(INDEX_3D_YZ, index + INDEX_Y) = GradNpT(i, INDEX_Z);
    result(INDEX_3D_YZ, index + INDEX_Z) = GradNpT(i, INDEX_Y);
    result(INDEX_3D_XZ, index + INDEX_X) = GradNpT(i, INDEX_Z);
    result(INDEX_3D_XZ, index + INDEX_Z) = GradNpT(i, INDEX_X);
    }

    return result;
}

double ThreeDimensionStressState::CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                                  double detJ,
                                                                  const Geometry<Node>&) const
{
    return rIntegrationPoint.Weight() * detJ;
}

unique_ptr<StressStatePolicy> ThreeDimensionStressState::Clone() const
{
    return std::make_unique<ThreeDimensionStressState>();
}

Vector ThreeDimensionStressState::CalculateGreenLagrangeStrain(const Matrix& rTotalDeformationGradient) const
{
    return Kratos::Vector();
}

}
