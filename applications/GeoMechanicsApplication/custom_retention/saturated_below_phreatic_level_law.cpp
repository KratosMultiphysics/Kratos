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
#include "geo_mechanics_application_variables.h"

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

    return 0;
}

} // Namespace Kratos
