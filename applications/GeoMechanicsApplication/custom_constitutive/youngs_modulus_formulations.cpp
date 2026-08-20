// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

#include <cmath>
#include <limits>
#include <type_traits>

namespace Kratos
{

namespace Formulations
{

YoungsModulusVariant InitializeFormulation(const std::string& rFormulation)
{
    if (rFormulation == Constant::Name) {
        return Constant{};
    } else if (rFormulation == Eur::Name) {
        return Eur{};
    } else {
        KRATOS_ERROR << "Unknown GEO_YOUNGS_MODULUS_FORMULATION: " << rFormulation;
    }
}

std::string GetYoungsModulusFormulation(const Properties& rProperties)
{
    return rProperties.Has(GEO_YOUNGS_MODULUS_FORMULATION) ? rProperties[GEO_YOUNGS_MODULUS_FORMULATION]
                                                           : Formulations::Constant::Name;
}

std::string GetYoungsModulusFormulation(const YoungsModulusVariant& rFormulation)
{
    return std::visit([](const auto& Formulation) -> std::string { return Formulation.Name; }, rFormulation);
}

double GetYoungsModulus(YoungsModulusVariant& rFormulation,
                        const Properties&     rProperties,
                        double                YoungsModulus,
                        Vector&               rStressVectorFinalized)
{
    return std::visit(
        [&rProperties, YoungsModulus, &rStressVectorFinalized]<typename TFormulation>(const TFormulation&) {
        if constexpr (std::is_same_v<std::decay_t<TFormulation>, Formulations::Constant>) {
            return YoungsModulus;

        } else if constexpr (std::is_same_v<std::decay_t<TFormulation>, Formulations::Eur>) {
            return Formulations::CalculateYoungsModulusForEur(rProperties, YoungsModulus, rStressVectorFinalized);
        }
    }, rFormulation);
}

double CalculateYoungsModulusForEur(const Properties& rProperties, double YoungsModulus, Vector& rStressVectorFinalized)
{
    constexpr auto epsilon = std::numeric_limits<double>::epsilon();

    const auto reference_pressure = rProperties[GEO_PRESSURE_REFERENCE];
    const auto exponent           = rProperties[GEO_STRESS_DEPENDENCY_EXPONENT];
    const auto eur_ref            = YoungsModulus;

    const auto friction_angle_rad = ConstitutiveLawUtilities::GetFrictionAngleInRadians(rProperties);
    const auto stress_shift =
        rProperties[GEO_COHESION] * std::cos(friction_angle_rad) / std::sin(friction_angle_rad);

    const auto base =
        (stress_shift - Formulations::CalculateMinorPrincipalEffectiveStress(rStressVectorFinalized)) /
        (stress_shift + reference_pressure);

    KRATOS_ERROR_IF_NOT(base > epsilon)
        << "Non-positive base for std::pow ("
        << base << "). Check GEO_COHESION, GEO_FRICTION_ANGLE, GEO_PRESSURE_REFERENCE and the finalized stress state.\n";

    return eur_ref * std::pow(base, exponent);
}

double CalculateMinorPrincipalEffectiveStress(Vector& rStressVectorFinalized)
{
    auto principal_stresses = Vector{};
    auto eigen_vectors      = Matrix{};
    StressStrainUtilities::CalculatePrincipalStresses(rStressVectorFinalized, principal_stresses, eigen_vectors);

    KRATOS_ERROR_IF(principal_stresses.size() < 3)
        << "Could not compute principal stresses from stress vector with size "
        << rStressVectorFinalized.size() << ". Expected at least 3 principal stresses, got "
        << principal_stresses.size() << "\n";

    return principal_stresses[2];
}

} // namespace Formulations
} // namespace Kratos
