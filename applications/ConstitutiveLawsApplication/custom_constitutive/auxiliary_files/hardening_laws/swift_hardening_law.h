// License: BSD License
// Kratos default license: kratos/license.txt

#pragma once

#include <cmath>

#include "includes/properties.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{

/** @brief Swift stress threshold and derivative with respect to accumulated plastic strain.
 * @details Stateless policy for GenericConstitutiveLawIntegratorPlasticity. It owns
 * neither history nor integration logic. Properties use the same stress units as E.
 */
class SwiftHardeningLaw
{
public:
    static void CalculateStressThreshold(
        const Properties& rProperties,
        const double AccumulatedPlasticStrain,
        double& rThreshold,
        double& rDerivative)
    {
        const double base = rProperties[SWIFT_INITIAL_STRAIN] + AccumulatedPlasticStrain;
        KRATOS_ERROR_IF(!std::isfinite(AccumulatedPlasticStrain) || AccumulatedPlasticStrain < 0.0 ||
            !std::isfinite(base) || base <= 0.0) << "Swift: invalid accumulated plastic strain or power base";
        const double exponent = rProperties[SWIFT_HARDENING_EXPONENT];
        rThreshold = rProperties[SWIFT_COEFFICIENT] * std::pow(base, exponent);
        rDerivative = exponent * (rThreshold / base);
        KRATOS_ERROR_IF(!std::isfinite(rThreshold) || rThreshold <= 0.0 ||
            !std::isfinite(rDerivative) || rDerivative <= 0.0)
            << "Swift: flow stress or hardening modulus is not representable; p = " << AccumulatedPlasticStrain;
    }

    static int Check(const Properties& rProperties)
    {
        for (const auto* p_variable : {&SWIFT_COEFFICIENT, &SWIFT_INITIAL_STRAIN, &SWIFT_HARDENING_EXPONENT}) {
            KRATOS_ERROR_IF_NOT(rProperties.Has(*p_variable)) << "Swift: missing " << p_variable->Name();
            const double value = rProperties[*p_variable];
            KRATOS_ERROR_IF(!std::isfinite(value) || value <= 0.0)
                << "Swift: " << p_variable->Name() << " must be finite and positive; got " << value;
        }
        double threshold = 0.0;
        double derivative = 0.0;
        CalculateStressThreshold(rProperties, 0.0, threshold, derivative);
        return 0;
    }
};

} // namespace Kratos
