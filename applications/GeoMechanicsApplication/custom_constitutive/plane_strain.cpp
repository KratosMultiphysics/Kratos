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
#include "includes/constitutive_law.h"

namespace Kratos
{

Matrix PlaneStrain::CalculateElasticConstitutiveTensor(const Properties& rProperties) const
{
    const double nu = false ? ConstitutiveLawUtilities::GetUndrainedPoissonsRatio(rProperties)
                            : rProperties[POISSON_RATIO];
    const double E  = false ? ConstitutiveLawUtilities::GetUndrainedYoungsModulus(rProperties, nu)
                            : rProperties[YOUNG_MODULUS];

    return ConstitutiveLawUtilities::MakeContinuumConstitutiveTensor(
        E, nu, PlaneStrain::GetStrainSize(), PlaneStrain::GetNumberOfNormalComponents());
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
