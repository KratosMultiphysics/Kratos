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
#include "custom_utilities/constitutive_law_utilities.h"
#include "geo_mechanics_application_constants.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

Matrix ThreeDimensional::CalculateElasticConstitutiveTensor(const Properties& rProperties) const
{
    constexpr auto undrained = false;

    const double nu = undrained ? ConstitutiveLawUtilities::GetUndrainedPoissonsRatio(rProperties)
                            : rProperties[POISSON_RATIO];
    const double E  = undrained ? ConstitutiveLawUtilities::GetUndrainedYoungsModulus(rProperties, nu)
                            : rProperties[YOUNG_MODULUS];

    return ConstitutiveLawUtilities::MakeContinuumConstitutiveTensor(
        E, nu, ThreeDimensional::GetStrainSize(), ThreeDimensional::GetNumberOfNormalComponents());
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
