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
#include "three_dimensional_stress_state.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "includes/serializer.h"

namespace Kratos
{

Matrix ThreeDimensionalStressState::CalculateBMatrix(const Matrix& rDN_DX,
                                                     const Vector&,
                                                     const Geometry<Node>& rGeometry) const
{
    const auto dimension       = rGeometry.WorkingSpaceDimension();
    const auto number_of_nodes = rGeometry.size();
    Matrix     result          = ZeroMatrix(VOIGT_SIZE_3D, dimension * number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        const auto offset = dimension * i;

        result(INDEX_3D_XX, offset + INDEX_X) = rDN_DX(i, INDEX_X);
        result(INDEX_3D_YY, offset + INDEX_Y) = rDN_DX(i, INDEX_Y);
        result(INDEX_3D_ZZ, offset + INDEX_Z) = rDN_DX(i, INDEX_Z);
        result(INDEX_3D_XY, offset + INDEX_X) = rDN_DX(i, INDEX_Y);
        result(INDEX_3D_XY, offset + INDEX_Y) = rDN_DX(i, INDEX_X);
        result(INDEX_3D_YZ, offset + INDEX_Y) = rDN_DX(i, INDEX_Z);
        result(INDEX_3D_YZ, offset + INDEX_Z) = rDN_DX(i, INDEX_Y);
        result(INDEX_3D_XZ, offset + INDEX_X) = rDN_DX(i, INDEX_Z);
        result(INDEX_3D_XZ, offset + INDEX_Z) = rDN_DX(i, INDEX_X);
    }

    return result;
}

Vector ThreeDimensionalStressState::CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const
{
    return MathUtils<>::StrainTensorToVector(
        StressStrainUtilities::CalculateGreenLagrangeStrainTensor(rDeformationGradient));
}

std::unique_ptr<StressStatePolicy> ThreeDimensionalStressState::Clone() const
{
    return std::make_unique<ThreeDimensionalStressState>();
}

const Vector& ThreeDimensionalStressState::GetVoigtVector() const { return VoigtVector3D; }

SizeType ThreeDimensionalStressState::GetVoigtSize() const { return GetVoigtSize3D(); }

SizeType ThreeDimensionalStressState::GetStressTensorSize() const
{
    return GetStressTensorSize3D();
}

void ThreeDimensionalStressState::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void ThreeDimensionalStressState::load(Serializer&)
{
    // No data members to be loaded (yet)
}

} // namespace Kratos
