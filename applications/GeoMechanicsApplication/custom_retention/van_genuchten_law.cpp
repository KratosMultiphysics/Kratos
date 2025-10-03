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
#include "custom_utilities/check_utilities.h"
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
    KRATOS_TRY
    const auto  p                     = rParameters.GetFluidPressure();
    const auto& r_material_properties = rParameters.GetMaterialProperties();

    if (p > 0.0) {
        const auto sat_max = r_material_properties[SATURATED_SATURATION];
        const auto sat_min = r_material_properties[RESIDUAL_SATURATION];
        const auto pb      = r_material_properties[VAN_GENUCHTEN_AIR_ENTRY_PRESSURE];
        const auto gn      = r_material_properties[VAN_GENUCHTEN_GN];
        const auto gc      = (1.0 - gn) / gn;

        return sat_min + (sat_max - sat_min) * std::pow(1.0 + std::pow(p / pb, gn), gc);
    } else {
        return r_material_properties[SATURATED_SATURATION];
    }

    KRATOS_CATCH("")
}

double VanGenuchtenLaw::CalculateEffectiveSaturation(Parameters& rParameters) const
{
    KRATOS_TRY

    const auto& r_material_properties = rParameters.GetMaterialProperties();
    const auto  sat_max               = r_material_properties[SATURATED_SATURATION];
    const auto  sat_min               = r_material_properties[RESIDUAL_SATURATION];

    return (CalculateSaturation(rParameters) - sat_min) / (sat_max - sat_min);

    KRATOS_CATCH("")
}

double VanGenuchtenLaw::CalculateDerivativeOfSaturation(Parameters& rParameters) const
{
    KRATOS_TRY
    const auto p = rParameters.GetFluidPressure();

    if (p > 0.0) {
        const auto& r_material_properties = rParameters.GetMaterialProperties();
        const auto  sat_max               = r_material_properties[SATURATED_SATURATION];
        const auto  sat_min               = r_material_properties[RESIDUAL_SATURATION];
        const auto  pb                    = r_material_properties[VAN_GENUCHTEN_AIR_ENTRY_PRESSURE];
        const auto  gn                    = r_material_properties[VAN_GENUCHTEN_GN];
        const auto  gc                    = (1.0 - gn) / gn;

        return (sat_max - sat_min) * gc * std::pow((1.0 + std::pow(p / pb, gn)), gc - 1.0) * gn *
               std::pow(pb, -gn) * std::pow(p, gn - 1.0);
    } else {
        return 0.0;
    }

    KRATOS_CATCH("")
}

double VanGenuchtenLaw::CalculateRelativePermeability(Parameters& rParameters) const
{
    KRATOS_TRY

    const auto eff_sat = CalculateEffectiveSaturation(rParameters);

    const auto& r_material_properties = rParameters.GetMaterialProperties();
    const auto  gl                    = r_material_properties[VAN_GENUCHTEN_GL];
    const auto  gn                    = r_material_properties[VAN_GENUCHTEN_GN];

    const auto rel_perm =
        std::pow(eff_sat, gl) *
        std::pow(1.0 - std::pow(1.0 - std::pow(eff_sat, gn / (gn - 1.0)), (gn - 1.0) / gn), 2);

    return std::max(rel_perm, r_material_properties[MINIMUM_RELATIVE_PERMEABILITY]);

    KRATOS_CATCH("")
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
    check_properties.Check(VAN_GENUCHTEN_GL);

    return 0;
}

std::string VanGenuchtenLaw::Info() const { return "VanGenuchtenLaw"s; }

} // namespace Kratos