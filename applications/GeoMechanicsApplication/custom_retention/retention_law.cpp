// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

/* Project includes */
#include "custom_retention/retention_law.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

double& RetentionLaw::CalculateValue(Parameters& rParameters, const Variable<double>& rThisVariable, double& rValue) const
{
    if (rThisVariable == DEGREE_OF_SATURATION) {
        rValue = CalculateSaturation(rParameters);
    } else if (rThisVariable == EFFECTIVE_SATURATION) {
        rValue = CalculateEffectiveSaturation(rParameters);
    } else if (rThisVariable == BISHOP_COEFFICIENT) {
        rValue = CalculateBishopCoefficient(rParameters);
    } else if (rThisVariable == DERIVATIVE_OF_SATURATION) {
        rValue = CalculateDerivativeOfSaturation(rParameters);
    } else if (rThisVariable == RELATIVE_PERMEABILITY) {
        rValue = CalculateRelativePermeability(rParameters);
    }

    return rValue;
}

void RetentionLaw::save(Serializer& rSerializer) const
{
    // there is no member variables to be saved
}

void RetentionLaw::load(Serializer& rSerializer)
{
    // there is no member variables to be loaded
}

} // namespace Kratos