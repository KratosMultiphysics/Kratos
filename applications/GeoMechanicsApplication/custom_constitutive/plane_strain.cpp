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

#include "plane_strain.h"

#include "geo_mechanics_application_constants.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

Matrix PlaneStrain::CalculateElasticMatrix(const Properties& rProperties) const
{
    const auto youngs_modulus = rProperties[YOUNG_MODULUS];
    const auto poissons_ratio = rProperties[POISSON_RATIO];
    const auto c0 = youngs_modulus / ((1.0 + poissons_ratio) * (1.0 - 2.0 * poissons_ratio));
    const auto c1 = (1.0 - poissons_ratio) * c0;
    const auto c2 = poissons_ratio * c0;
    const auto c3 = (0.5 - poissons_ratio) * c0;

    Matrix result = ZeroMatrix(4, 4);

    result(INDEX_2D_PLANE_STRAIN_XX, INDEX_2D_PLANE_STRAIN_XX) = c1;
    result(INDEX_2D_PLANE_STRAIN_XX, INDEX_2D_PLANE_STRAIN_YY) = c2;
    result(INDEX_2D_PLANE_STRAIN_XX, INDEX_2D_PLANE_STRAIN_ZZ) = c2;

    result(INDEX_2D_PLANE_STRAIN_YY, INDEX_2D_PLANE_STRAIN_XX) = c2;
    result(INDEX_2D_PLANE_STRAIN_YY, INDEX_2D_PLANE_STRAIN_YY) = c1;
    result(INDEX_2D_PLANE_STRAIN_YY, INDEX_2D_PLANE_STRAIN_ZZ) = c2;

    result(INDEX_2D_PLANE_STRAIN_ZZ, INDEX_2D_PLANE_STRAIN_XX) = c2;
    result(INDEX_2D_PLANE_STRAIN_ZZ, INDEX_2D_PLANE_STRAIN_YY) = c2;
    result(INDEX_2D_PLANE_STRAIN_ZZ, INDEX_2D_PLANE_STRAIN_ZZ) = c1;

    result(INDEX_2D_PLANE_STRAIN_XY, INDEX_2D_PLANE_STRAIN_XY) = c3;

    return result;
}

std::unique_ptr<ConstitutiveLawDimension> PlaneStrain::Clone() const
{
    return std::make_unique<PlaneStrain>();
}

std::size_t PlaneStrain::GetStrainSize() const { return VOIGT_SIZE_2D_PLANE_STRAIN; }

std::size_t PlaneStrain::GetDimension() const { return N_DIM_2D; }

std::size_t PlaneStrain::GetNumberOfNormalComponents() const { return 3; }

Flags PlaneStrain::GetSpatialType() const { return ConstitutiveLaw::PLANE_STRAIN_LAW; }

void PlaneStrain::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void PlaneStrain::load(Serializer&)
{
    // No data members to be saved (yet)
}

} // namespace Kratos
