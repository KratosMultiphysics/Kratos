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
//                   Gennady Markelov
//

#include "three_dimensional.h"

#include "geo_mechanics_application_constants.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

Matrix ThreeDimensional::CalculateElasticMatrix(double YoungsModulus, double PoissonsRatio) const
{
    const auto c0 = YoungsModulus / ((1.0 + PoissonsRatio) * (1.0 - 2.0 * PoissonsRatio));
    const auto c1 = (1.0 - PoissonsRatio) * c0;
    const auto c2 = PoissonsRatio * c0;
    const auto c3 = (0.5 - PoissonsRatio) * c0;

    Matrix result = ZeroMatrix(GetStrainSize(), GetStrainSize());

    result(INDEX_3D_XX, INDEX_3D_XX) = c1;
    result(INDEX_3D_XX, INDEX_3D_YY) = c2;
    result(INDEX_3D_XX, INDEX_3D_ZZ) = c2;

    result(INDEX_3D_YY, INDEX_3D_XX) = c2;
    result(INDEX_3D_YY, INDEX_3D_YY) = c1;
    result(INDEX_3D_YY, INDEX_3D_ZZ) = c2;

    result(INDEX_3D_ZZ, INDEX_3D_XX) = c2;
    result(INDEX_3D_ZZ, INDEX_3D_YY) = c2;
    result(INDEX_3D_ZZ, INDEX_3D_ZZ) = c1;

    result(INDEX_3D_XY, INDEX_3D_XY) = c3;
    result(INDEX_3D_YZ, INDEX_3D_YZ) = c3;
    result(INDEX_3D_XZ, INDEX_3D_XZ) = c3;

    return result;
}

std::unique_ptr<ConstitutiveLawDimension> ThreeDimensional::Clone() const
{
    return std::make_unique<ThreeDimensional>();
}

std::size_t ThreeDimensional::GetStrainSize() const { return VOIGT_SIZE_3D; }

std::size_t ThreeDimensional::GetDimension() const { return N_DIM_3D; }

std::size_t ThreeDimensional::GetNumberOfNormalComponents() const { return 3; }

Flags ThreeDimensional::GetSpatialType() const { return ConstitutiveLaw::THREE_DIMENSIONAL_LAW; }

void ThreeDimensional::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void ThreeDimensional::load(Serializer&)
{
    // No data members to be saved (yet)
}

} // namespace Kratos
