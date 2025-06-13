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
#include "geo_mechanics_application_variables.h"

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
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SATURATED_SATURATION))
        << "SATURATED_SATURATION is not available in the parameters of material "
        << rMaterialProperties.Id() << "." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] < 0.0 || rMaterialProperties[SATURATED_SATURATION] > 1.0)
        << "SATURATED_SATURATION (" << rMaterialProperties[SATURATED_SATURATION]
        << ") must be in the range [0.0, 1.0] for material " << rMaterialProperties.Id() << "."
        << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(RESIDUAL_SATURATION))
        << "RESIDUAL_SATURATION is not available in the parameters of material "
        << rMaterialProperties.Id() << "." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[RESIDUAL_SATURATION] < 0.0 || rMaterialProperties[RESIDUAL_SATURATION] >= rMaterialProperties[SATURATED_SATURATION])
        << "RESIDUAL_SATURATION (" << rMaterialProperties[RESIDUAL_SATURATION]
        << ") must be in the range [0.0, " << rMaterialProperties[SATURATED_SATURATION]
        << "> for material " << rMaterialProperties.Id() << "." << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(MINIMUM_RELATIVE_PERMEABILITY))
        << "MINIMUM_RELATIVE_PERMEABILITY is not available in the parameters of material "
        << rMaterialProperties.Id() << "." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[MINIMUM_RELATIVE_PERMEABILITY] < 0.0 || rMaterialProperties[MINIMUM_RELATIVE_PERMEABILITY] > 1.0)
        << "MINIMUM_RELATIVE_PERMEABILITY (" << rMaterialProperties[MINIMUM_RELATIVE_PERMEABILITY]
        << ") must be in the range [0.0, 1.0] for material " << rMaterialProperties.Id() << "."
        << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(VAN_GENUCHTEN_AIR_ENTRY_PRESSURE))
        << "VAN_GENUCHTEN_AIR_ENTRY_PRESSURE is not available in the parameters of material "
        << rMaterialProperties.Id() << "." << std::endl;
    KRATOS_ERROR_IF_NOT((rMaterialProperties[VAN_GENUCHTEN_AIR_ENTRY_PRESSURE] > 0.0))
        << "VAN_GENUCHTEN_AIR_ENTRY_PRESSURE ("
        << rMaterialProperties[VAN_GENUCHTEN_AIR_ENTRY_PRESSURE] << ") must be greater than 0 "
        << "for material " << rMaterialProperties.Id() << "." << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(VAN_GENUCHTEN_GN))
        << "VAN_GENUCHTEN_GN is not available in the parameters of material "
        << rMaterialProperties.Id() << "." << std::endl;
    KRATOS_ERROR_IF_NOT((rMaterialProperties[VAN_GENUCHTEN_GN] > 0.0))
        << "VAN_GENUCHTEN_GN (" << rMaterialProperties[VAN_GENUCHTEN_GN] << ") must be greater than 0 "
        << "for material " << rMaterialProperties.Id() << "." << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(VAN_GENUCHTEN_GL))
        << "VAN_GENUCHTEN_GL is not available in the parameters of material "
        << rMaterialProperties.Id() << "." << std::endl;
    KRATOS_ERROR_IF_NOT((rMaterialProperties[VAN_GENUCHTEN_GL] >= 0.0))
        << "VAN_GENUCHTEN_GL (" << rMaterialProperties[VAN_GENUCHTEN_GL]
        << ") must be greater than or equal to 0 for material " << rMaterialProperties.Id() << "."
        << std::endl;

    return 0;
}

} // namespace Kratos