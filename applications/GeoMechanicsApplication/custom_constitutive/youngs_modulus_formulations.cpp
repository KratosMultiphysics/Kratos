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

#include "custom_constitutive/youngs_modulus_formulations.h"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

#include <cmath>
#include <limits>
#include <type_traits>

namespace Kratos
{
double GeoYoungsModulusFormulations::Constant::operator()(const Properties&,
                                                          double ReferenceYoungsModulus,
                                                          const Vector&) const
{
    return ReferenceYoungsModulus;
}

double GeoYoungsModulusFormulations::SchanzVermeer::operator()(const Properties& rProperties,
                                                               double ReferenceYoungsModulus,
                                                               const Vector& rStressVectorFinalized) const
{
    return CalculateYoungsModulusSchanzVermeer(rProperties, ReferenceYoungsModulus, rStressVectorFinalized);
}

void GeoYoungsModulusFormulations::CheckInputData(const Properties& rMaterialProperties)
{
    if (rMaterialProperties.Has(GEO_YOUNGS_MODULUS_FORMULATION) &&
        rMaterialProperties[GEO_YOUNGS_MODULUS_FORMULATION] == SchanzVermeer::Name) {
        const CheckProperties check_properties(rMaterialProperties, "parameters of material",
                                               CheckProperties::Bounds::AllExclusive);
        check_properties.Check(GEO_PRESSURE_REFERENCE);
        check_properties.Check(GEO_STRESS_DEPENDENCY_EXPONENT);
        check_properties.Check(GEO_COHESION);
        check_properties.Check(GEO_FRICTION_ANGLE);
    }
}

GeoYoungsModulusFormulations::YoungsModulusVariant GeoYoungsModulusFormulations::InitializeFormulation(const std::string& rFormulation)
{
    if (rFormulation == Constant::Name) return Constant{};
    if (rFormulation == SchanzVermeer::Name) return SchanzVermeer{};
    KRATOS_ERROR << "Unknown GEO_YOUNGS_MODULUS_FORMULATION: " << rFormulation;
}

std::string GeoYoungsModulusFormulations::GetYoungsModulusFormulation(const Properties& rProperties)
{
    return rProperties.Has(GEO_YOUNGS_MODULUS_FORMULATION) ? rProperties[GEO_YOUNGS_MODULUS_FORMULATION]
                                                           : Constant::Name;
}

std::string GeoYoungsModulusFormulations::GetYoungsModulusFormulation(const YoungsModulusVariant& rFormulation)
{
    return std::visit([](const auto& Formulation) -> std::string { return Formulation.Name; }, rFormulation);
}

double GeoYoungsModulusFormulations::GetYoungsModulus(YoungsModulusVariant& rFormulation,
                                                      const Properties&     rProperties,
                                                      double                ReferenceYoungsModulus,
                                                      const Vector&         rStressVectorFinalized)
{
    return std::visit([&](const auto& formulation_functor) {
        return formulation_functor(rProperties, ReferenceYoungsModulus, rStressVectorFinalized);
    }, rFormulation);
}

double GeoYoungsModulusFormulations::CalculateYoungsModulusSchanzVermeer(const Properties& rProperties,
                                                                         double ReferenceYoungsModulus,
                                                                         const Vector& rStressVectorFinalized)
{
    constexpr auto epsilon = std::numeric_limits<double>::epsilon();

    const auto reference_pressure = rProperties[GEO_PRESSURE_REFERENCE];
    const auto exponent           = rProperties[GEO_STRESS_DEPENDENCY_EXPONENT];

    const auto friction_angle_rad = ConstitutiveLawUtilities::GetFrictionAngleInRadians(rProperties);
    const auto stress_shift = rProperties[GEO_COHESION] / std::tan(friction_angle_rad);

    const auto base = (stress_shift - CalculateMinorPrincipalEffectiveStress(rStressVectorFinalized)) /
                      (stress_shift + reference_pressure);

    KRATOS_ERROR_IF_NOT(base > epsilon)
        << "Non-positive base for std::pow ("
        << base << "). Check GEO_COHESION, GEO_FRICTION_ANGLE, GEO_PRESSURE_REFERENCE and the finalized stress state.\n";

    return ReferenceYoungsModulus * std::pow(base, exponent);
}

double GeoYoungsModulusFormulations::CalculateMinorPrincipalEffectiveStress(const Vector& rStressVectorFinalized)
{
    auto principal_stresses = Vector{};
    auto eigen_vectors      = Matrix{};
    StressStrainUtilities::CalculatePrincipalStresses(rStressVectorFinalized, principal_stresses, eigen_vectors);

    KRATOS_ERROR_IF(principal_stresses.size() != 3)
        << "Could not compute principal stresses from stress vector with size "
        << rStressVectorFinalized.size() << ". Expected 3 principal stresses, got "
        << principal_stresses.size() << "." << "\n";

    return principal_stresses[2];
}

} // namespace Kratos
