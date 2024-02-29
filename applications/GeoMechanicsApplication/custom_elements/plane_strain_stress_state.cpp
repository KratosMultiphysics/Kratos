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
#include "custom_utilities/stress_strain_utilities.h"

namespace Kratos
{

Matrix PlaneStrainStressState::CalculateBMatrix(const Matrix& rGradNpT, const Vector&, const Geometry<Node>& rGeometry) const
{
    const auto dimension       = rGeometry.WorkingSpaceDimension();
    const auto number_of_nodes = rGeometry.size();
    Matrix     result = ZeroMatrix(VOIGT_SIZE_2D_AXISYMMETRIC, dimension * number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        const auto offset = dimension * i;

        result(INDEX_2D_PLANE_STRAIN_XX, offset + INDEX_X) = rGradNpT(i, INDEX_X);
        result(INDEX_2D_PLANE_STRAIN_YY, offset + INDEX_Y) = rGradNpT(i, INDEX_Y);
        result(INDEX_2D_PLANE_STRAIN_XY, offset + INDEX_X) = rGradNpT(i, INDEX_Y);
        result(INDEX_2D_PLANE_STRAIN_XY, offset + INDEX_Y) = rGradNpT(i, INDEX_X);
    }

    return result;
}

double PlaneStrainStressState::CalculateIntegrationCoefficient(const Geometry<Node>::IntegrationPointType& rIntegrationPoint,
                                                               double DetJ,
                                                               const Geometry<Node>&) const
{
    return rIntegrationPoint.Weight() * DetJ;
}

Vector PlaneStrainStressState::CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const
{
    return ConvertStrainTensorToVector(StressStrainUtilities::CalculateGreenLagrangeStrainTensor(rDeformationGradient));
}

std::unique_ptr<StressStatePolicy> PlaneStrainStressState::Clone() const
{
    return std::make_unique<PlaneStrainStressState>();
}

Vector PlaneStrainStressState::ConvertStrainTensorToVector(const Matrix& rStrainTensor)
{
    const auto strain_vector         = MathUtils<double>::StrainTensorToVector(rStrainTensor);
    Vector     result                = ZeroVector(VOIGT_SIZE_2D_PLANE_STRAIN);
    result[INDEX_2D_PLANE_STRAIN_XX] = strain_vector[0];
    result[INDEX_2D_PLANE_STRAIN_YY] = strain_vector[1];
    result[INDEX_2D_PLANE_STRAIN_ZZ] = 0.0;
    result[INDEX_2D_PLANE_STRAIN_XY] = strain_vector[2];
    return result;
}

} // namespace Kratos
