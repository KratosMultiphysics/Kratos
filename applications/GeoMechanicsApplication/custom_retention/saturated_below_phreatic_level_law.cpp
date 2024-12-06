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

// System includes
#include "custom_retention/saturated_below_phreatic_level_law.h"
#include <iostream>

namespace Kratos
{
SaturatedBelowPhreaticLevelLaw::SaturatedBelowPhreaticLevelLaw() : RetentionLaw() {}

SaturatedBelowPhreaticLevelLaw::SaturatedBelowPhreaticLevelLaw(const SaturatedBelowPhreaticLevelLaw& rOther)
    : RetentionLaw(rOther)
{
}

RetentionLaw::Pointer SaturatedBelowPhreaticLevelLaw::Clone() const
{
    return Kratos::make_shared<SaturatedBelowPhreaticLevelLaw>(*this);
}

SaturatedBelowPhreaticLevelLaw::~SaturatedBelowPhreaticLevelLaw() {}

double SaturatedBelowPhreaticLevelLaw::CalculateSaturation(Parameters& rParameters) const
{
    if (rParameters.GetFluidPressure() < 0.0) {
        return rParameters.GetMaterialProperties()[SATURATED_SATURATION];
    } else {
        return rParameters.GetMaterialProperties()[RESIDUAL_SATURATION];
    }
}

double SaturatedBelowPhreaticLevelLaw::CalculateEffectiveSaturation(Parameters& rParameters) const
{
    const auto& r_material_properties = rParameters.GetMaterialProperties();
    const auto& sat_max              = r_material_properties[SATURATED_SATURATION];
    const auto& sat_min              = r_material_properties[RESIDUAL_SATURATION];

    return (CalculateSaturation(rParameters) - sat_min) / (sat_max - sat_min);
}

double SaturatedBelowPhreaticLevelLaw::CalculateDerivativeOfSaturation(Parameters& rParameters) const
{
    return 0.0;
}

double SaturatedBelowPhreaticLevelLaw::CalculateRelativePermeability(Parameters& rParameters) const
{
    if (rParameters.GetFluidPressure() < 0.0) {
        return 1.0;
    } else {
        return rParameters.GetMaterialProperties()[MINIMUM_RELATIVE_PERMEABILITY];
    }
}

double SaturatedBelowPhreaticLevelLaw::CalculateBishopCoefficient(Parameters& rParameters) const
{
    return CalculateEffectiveSaturation(rParameters);
}

double& SaturatedBelowPhreaticLevelLaw::CalculateValue(Parameters&             rParameterValues,
                                                       const Variable<double>& rThisVariable,
                                                       double&                 rValue)
{
    if (rThisVariable == DEGREE_OF_SATURATION) {
        rValue = this->CalculateSaturation(rParameterValues);
    } else if (rThisVariable == EFFECTIVE_SATURATION) {
        rValue = this->CalculateEffectiveSaturation(rParameterValues);
    } else if (rThisVariable == BISHOP_COEFFICIENT) {
        rValue = this->CalculateBishopCoefficient(rParameterValues);
    } else if (rThisVariable == DERIVATIVE_OF_SATURATION) {
        rValue = this->CalculateDerivativeOfSaturation(rParameterValues);
    } else if (rThisVariable == RELATIVE_PERMEABILITY) {
        rValue = this->CalculateRelativePermeability(rParameterValues);
    }
    return rValue;
}

void SaturatedBelowPhreaticLevelLaw::InitializeMaterial(const Properties&, const GeometryType&, const Vector&)
{
    // nothing is needed
}

void SaturatedBelowPhreaticLevelLaw::Initialize(Parameters&)
{
    // nothing is needed
}

void SaturatedBelowPhreaticLevelLaw::InitializeSolutionStep(Parameters&)
{
    // nothing is needed
}

void SaturatedBelowPhreaticLevelLaw::Finalize(Parameters&)
{
    // nothing is needed
}

void SaturatedBelowPhreaticLevelLaw::FinalizeSolutionStep(Parameters&)
{
    // nothing is needed
}

int SaturatedBelowPhreaticLevelLaw::Check(const Properties& rMaterialProperties, const ProcessInfo&)
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SATURATED_SATURATION))
        << "SATURATED_SATURATION is not available in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] < 0.0)
        << "SATURATED_SATURATION cannot be less than 0 " << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(RESIDUAL_SATURATION))
        << "RESIDUAL_SATURATION is not available in material parameters" << std::endl;
    KRATOS_DEBUG_ERROR_IF_NOT(rMaterialProperties[RESIDUAL_SATURATION] > 0.0)
        << "RESIDUAL_SATURATION must be greater than 0 " << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[RESIDUAL_SATURATION] > 1.0)
        << "RESIDUAL_SATURATION cannot be greater than 1.0 " << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] < rMaterialProperties[RESIDUAL_SATURATION])
        << "RESIDUAL_SATURATION cannot be greater than SATURATED_SATURATION " << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(MINIMUM_RELATIVE_PERMEABILITY))
        << "MINIMUM_RELATIVE_PERMEABILITY is not available in material parameters" << std::endl;
    KRATOS_ERROR_IF_NOT((rMaterialProperties[MINIMUM_RELATIVE_PERMEABILITY] > 0.0))
        << "MINIMUM_RELATIVE_PERMEABILITY must be greater than 0 " << std::endl;

    return 0;
}

} // Namespace Kratos
