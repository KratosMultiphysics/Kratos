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
#include <iostream>

#include "custom_retention/saturated_law.h"

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

double& SaturatedLaw::CalculateValue(RetentionLaw::Parameters& rParameterValues,
                                     const Variable<double>&   rThisVariable,
                                     double&                   rValue)
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

void SaturatedLaw::InitializeMaterial(const Properties&   rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector&       rShapeFunctionsValues)
{
    // nothing is needed
}

void SaturatedLaw::Initialize(Parameters& rParameters)
{
    // nothing is needed
}

void SaturatedLaw::InitializeSolutionStep(Parameters& rParameters)
{
    // nothing is needed
}

void SaturatedLaw::Finalize(Parameters& rParameters)
{
    // nothing is needed
}

void SaturatedLaw::FinalizeSolutionStep(Parameters& rParameters)
{
    // nothing is needed
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