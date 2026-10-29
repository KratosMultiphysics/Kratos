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

#include "custom_retention/van_genuchten_law.h"
#include "custom_utilities/check_utilities.hpp"
#include "geo_mechanics_application_variables.h"

#include <string>

using namespace std::string_literals;

namespace Kratos
{

RetentionLaw::Pointer VanGenuchtenLaw::Clone() const
{
    return Kratos::make_shared<VanGenuchtenLaw>(*this);
}

double VanGenuchtenLaw::CalculateSaturation(Parameters& rParameters) const
{
    const auto& r_material_properties = rParameters.GetMaterialProperties();

    const auto s_s = r_material_properties[SATURATED_SATURATION];
    const auto s_r = r_material_properties[RESIDUAL_SATURATION];

    return s_r + CalculateEffectiveSaturation(rParameters) * (s_s - s_r);
}

double VanGenuchtenLaw::CalculateEffectiveSaturation(Parameters& rParameters) const
{
    if (const auto p = rParameters.GetFluidPressure(); p > 0.0) {
        const auto& r_material_properties = rParameters.GetMaterialProperties();
        const auto  p_b                   = r_material_properties[VAN_GENUCHTEN_AIR_ENTRY_PRESSURE];
        const auto  n                     = r_material_properties[VAN_GENUCHTEN_GN];
        const auto  m                     = (n - 1.0) / n;

        return std::pow(1.0 + std::pow(p / p_b, n), -m);
    }
    return 1.0;
}

double VanGenuchtenLaw::CalculateDerivativeOfSaturation(Parameters& rParameters) const
{
    if (const auto p = rParameters.GetFluidPressure(); p > 0.0) {
        const auto& r_material_properties = rParameters.GetMaterialProperties();
        const auto  s_s                   = r_material_properties[SATURATED_SATURATION];
        const auto  s_r                   = r_material_properties[RESIDUAL_SATURATION];
        const auto  p_b                   = r_material_properties[VAN_GENUCHTEN_AIR_ENTRY_PRESSURE];
        const auto  n                     = r_material_properties[VAN_GENUCHTEN_GN];
        const auto  m                     = (n - 1.0) / n;

        return (s_s - s_r) * (-m) * n * std::pow(p_b, -n) * std::pow(p, n - 1.0) /
               std::pow((1.0 + std::pow(p / p_b, n)), m + 1.0);
    }
    return 0.0;
}

double VanGenuchtenLaw::CalculateRelativePermeability(Parameters& rParameters) const
{
    const auto eff_sat = CalculateEffectiveSaturation(rParameters);

    const auto& r_material_properties = rParameters.GetMaterialProperties();
    const auto  l                     = r_material_properties[VAN_GENUCHTEN_GL];
    const auto  n                     = r_material_properties[VAN_GENUCHTEN_GN];
    const auto  m                     = (n - 1.0) / n;

    const auto rel_perm =
        std::pow(eff_sat, l) * std::pow(1.0 - std::pow(1.0 - std::pow(eff_sat, 1.0 / m), m), 2);

    return std::max(rel_perm, r_material_properties[MINIMUM_RELATIVE_PERMEABILITY]);
}

double VanGenuchtenLaw::CalculateBishopCoefficient(Parameters& rParameters) const
{
    return CalculateEffectiveSaturation(rParameters);
}

int VanGenuchtenLaw::Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo)
{
    using enum CheckProperties::Bounds;
    const CheckProperties check_properties(rMaterialProperties, "parameters of material", AllInclusive);
    constexpr auto max_value = 1.0;
    check_properties.Check(SATURATED_SATURATION, max_value);
    check_properties.SingleUseBounds(InclusiveLowerAndExclusiveUpper)
        .Check(RESIDUAL_SATURATION, rMaterialProperties[SATURATED_SATURATION]);
    check_properties.Check(MINIMUM_RELATIVE_PERMEABILITY, max_value);
    check_properties.SingleUseBounds(AllExclusive).Check(VAN_GENUCHTEN_AIR_ENTRY_PRESSURE);
    check_properties.SingleUseBounds(AllExclusive).Check(VAN_GENUCHTEN_GN);
    check_properties.CheckAvailability(VAN_GENUCHTEN_GL);

    return 0;
}

std::string VanGenuchtenLaw::Info() const { return "VanGenuchtenLaw"s; }

} // namespace Kratos
