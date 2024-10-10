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

Matrix PlaneStrain::FillConstitutiveMatrix(double c1, double c2, double c3) const
{
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

Flags PlaneStrain::GetSpatialType() const { return ConstitutiveLaw::PLANE_STRAIN_LAW; }

} // namespace Kratos
