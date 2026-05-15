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
#include "includes/serializer.h"

namespace Kratos
{
using enum indexDOF3D;
using enum indexStress2DPlaneStrain;

Matrix PlaneStrainStressState::CalculateBMatrix(const Matrix& rDN_DX, const Vector&, const Geometry<Node>& rGeometry) const
{
    const auto dimension       = rGeometry.WorkingSpaceDimension();
    const auto number_of_nodes = rGeometry.size();
    Matrix     result = ZeroMatrix(VOIGT_SIZE_2D_PLANE_STRAIN, dimension * number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        const auto offset = dimension * i;

        result(static_cast<std::size_t>(INDEX_2D_PLANE_STRAIN_XX),
               offset + static_cast<std::size_t>(INDEX_X)) = rDN_DX(i, static_cast<std::size_t>(INDEX_X));
        result(static_cast<std::size_t>(INDEX_2D_PLANE_STRAIN_YY),
               offset + static_cast<std::size_t>(INDEX_Y)) = rDN_DX(i, static_cast<std::size_t>(INDEX_Y));
        result(static_cast<std::size_t>(INDEX_2D_PLANE_STRAIN_XY),
               offset + static_cast<std::size_t>(INDEX_X)) = rDN_DX(i, static_cast<std::size_t>(INDEX_Y));
        result(static_cast<std::size_t>(INDEX_2D_PLANE_STRAIN_XY),
               offset + static_cast<std::size_t>(INDEX_Y)) = rDN_DX(i, static_cast<std::size_t>(INDEX_X));
    }

    return result;
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
    const auto strain_vector = MathUtils<double>::StrainTensorToVector(rStrainTensor);
    Vector     result        = ZeroVector(VOIGT_SIZE_2D_PLANE_STRAIN);
    result[static_cast<std::size_t>(INDEX_2D_PLANE_STRAIN_XX)] = strain_vector[0];
    result[static_cast<std::size_t>(INDEX_2D_PLANE_STRAIN_YY)] = strain_vector[1];
    result[static_cast<std::size_t>(INDEX_2D_PLANE_STRAIN_ZZ)] = 0.0;
    result[static_cast<std::size_t>(INDEX_2D_PLANE_STRAIN_XY)] = strain_vector[2];
    return result;
}

const Vector& PlaneStrainStressState::GetVoigtVector() const { return VoigtVector2D; }

SizeType PlaneStrainStressState::GetVoigtSize() const { return GetVoigtSize2D(); }

SizeType PlaneStrainStressState::GetStressTensorSize() const { return GetStressTensorSize2D(); }

void PlaneStrainStressState::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void PlaneStrainStressState::load(Serializer&)
{
    // No data members to be loaded (yet)
}

} // namespace Kratos
