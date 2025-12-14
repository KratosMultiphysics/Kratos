// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Vahid Galavi
//

#include "custom_retention/saturated_below_phreatic_level_law.h"
#include "custom_utilities/check_utilities.h"
#include "geo_mechanics_application_variables.h"

#include <string>

using namespace std::string_literals;

namespace Kratos
{
RetentionLaw::Pointer SaturatedBelowPhreaticLevelLaw::Clone() const
{
    return Kratos::make_shared<SaturatedBelowPhreaticLevelLaw>(*this);
}

double SaturatedBelowPhreaticLevelLaw::CalculateSaturation(Parameters& rParameters) const
{
    return (rParameters.GetFluidPressure() < 0.0)
               ? rParameters.GetMaterialProperties()[SATURATED_SATURATION]
               : rParameters.GetMaterialProperties()[RESIDUAL_SATURATION];
}

double SaturatedBelowPhreaticLevelLaw::CalculateEffectiveSaturation(Parameters& rParameters) const
{
    const auto& r_material_properties = rParameters.GetMaterialProperties();
    const auto  sat_max               = r_material_properties[SATURATED_SATURATION];
    const auto  sat_min               = r_material_properties[RESIDUAL_SATURATION];

    return (CalculateSaturation(rParameters) - sat_min) / (sat_max - sat_min);
}

double SaturatedBelowPhreaticLevelLaw::CalculateDerivativeOfSaturation(Parameters&) const
{
    return 0.0;
}

double SaturatedBelowPhreaticLevelLaw::CalculateRelativePermeability(Parameters& rParameters) const
{
    return rParameters.GetFluidPressure() < 0.0
               ? 1.0
               : rParameters.GetMaterialProperties()[MINIMUM_RELATIVE_PERMEABILITY];
}

double SaturatedBelowPhreaticLevelLaw::CalculateBishopCoefficient(Parameters& rParameters) const
{
    return CalculateEffectiveSaturation(rParameters);
}

int SaturatedBelowPhreaticLevelLaw::Check(const Properties& rMaterialProperties, const ProcessInfo&)
{
    const CheckProperties check_properties(rMaterialProperties, "parameters of material",
                                           CheckProperties::Bounds::AllInclusive);
    constexpr auto        max_value = 1.0;
    check_properties.Check(SATURATED_SATURATION, max_value);
    check_properties.SingleUseBounds(CheckProperties::Bounds::InclusiveLowerAndExclusiveUpper)
        .Check(RESIDUAL_SATURATION, rMaterialProperties[SATURATED_SATURATION]);
    check_properties.Check(MINIMUM_RELATIVE_PERMEABILITY, max_value);

    return 0;
}

std::string SaturatedBelowPhreaticLevelLaw::Info() const
{
    return "SaturatedBelowPhreaticLevelLaw"s;
}
} // Namespace Kratos
