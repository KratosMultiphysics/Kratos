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

#include "custom_retention/saturated_law.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

RetentionLaw::Pointer SaturatedLaw::Clone() const
{
    return Kratos::make_shared<SaturatedLaw>(*this);
}

double SaturatedLaw::CalculateSaturation(Parameters& rParameters) const
{
    const Properties& rMaterialProperties = rParameters.GetMaterialProperties();
    return rMaterialProperties.Has(SATURATED_SATURATION) ? rMaterialProperties[SATURATED_SATURATION] : 1.0;
}

double SaturatedLaw::CalculateEffectiveSaturation(Parameters& rParameters) const { return 1.0; }

double SaturatedLaw::CalculateDerivativeOfSaturation(Parameters& rParameters) const { return 0.0; }

double SaturatedLaw::CalculateRelativePermeability(Parameters& rParameters) const { return 1.0; }

double SaturatedLaw::CalculateBishopCoefficient(Parameters& rParameters) const
{
    return CalculateEffectiveSaturation(rParameters);
}

int SaturatedLaw::Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo)
{
    if (rMaterialProperties.Has(SATURATED_SATURATION)) {
        KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] < 0.0 || rMaterialProperties[SATURATED_SATURATION] > 1.0)
            << "SATURATED_SATURATION (" << rMaterialProperties[SATURATED_SATURATION]
            << ") must be in the range [0.0, 1.0] for material " << rMaterialProperties.Id() << "."
            << std::endl;
    }

    return 0;
}

} // namespace Kratos