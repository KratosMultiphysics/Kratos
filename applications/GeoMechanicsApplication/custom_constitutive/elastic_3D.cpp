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

#include "elastic_3D.h"

#include "geo_mechanics_application_constants.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

Matrix Elastic3D::FillConstitutiveMatrix(double c1, double c2, double c3) const
{
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

std::unique_ptr<ConstitutiveLawDimension> Elastic3D::Clone() const
{
    return std::make_unique<Elastic3D>();
}

std::size_t Elastic3D::GetStrainSize() const { return VOIGT_SIZE_3D; }

std::size_t Elastic3D::GetDimension() const { return N_DIM_3D; }

Flags Elastic3D::GetSpatialType() const { return ConstitutiveLaw::THREE_DIMENSIONAL_LAW; }

} // namespace Kratos
