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

#include "custom_utilities/constitutive_law_utilities.h"
#include "geo_mechanics_application_constants.h"
#include "geo_mechanics_application_variables.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

Matrix PlaneStrain::CalculateElasticMatrix(const Properties& rProperties) const
{
    double E;
    double nu;
    auto   drainage_type = rProperties.Has(GEO_DRAINAGE_TYPE)
                               ? static_cast<DrainageType>(rProperties[GEO_DRAINAGE_TYPE])
                               : DrainageType::FULLY_COUPLED;
    if (drainage_type == DrainageType::UNDRAINED) {
        nu = ConstitutiveLawUtilities::GetUndrainedPoissonsRatio(rProperties);
        E  = ConstitutiveLawUtilities::GetUndrainedYoungsModulus(rProperties, nu);
    } else {
        E  = rProperties[YOUNG_MODULUS];
        nu = rProperties[POISSON_RATIO];
    }

    const auto c0            = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const auto c1            = (1.0 - nu) * c0;
    const auto c2            = nu * c0;
    const auto shear_modulus = E / (2.0 * (1.0 + nu));

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

    result(INDEX_2D_PLANE_STRAIN_XY, INDEX_2D_PLANE_STRAIN_XY) = shear_modulus;

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
