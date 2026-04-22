// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//                   Wijtze Pieter Kikstra
//

#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/string_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "utilities/math_utils.h"
#include <algorithm>

using namespace std::string_literals;

namespace
{

using namespace Kratos;

double GetValueOfUMatParameter(const Properties& rProperties, const Variable<int>& rIndexVariable)
{
    KRATOS_ERROR_IF_NOT(rProperties.Has(UMAT_PARAMETERS))
        << "There is no UMAT_PARAMETERS for material " << rProperties.Id() << "." << std::endl;

    KRATOS_ERROR_IF_NOT(rProperties.Has(rIndexVariable))
        << "There is no " << rIndexVariable.Name() << " for material " << rProperties.Id() << "."
        << std::endl;

    const auto index = rProperties[rIndexVariable]; // 1-based index
    KRATOS_DEBUG_ERROR_IF(index < 1 ||
                          static_cast<std::size_t>(index) > rProperties[UMAT_PARAMETERS].size())
        << "Got out-of-bounds " << rIndexVariable.Name() << " (material ID: " << rProperties.Id()
        << "): " << index << " is not in range [1, " << rProperties[UMAT_PARAMETERS].size() << "].\n";

    return rProperties[UMAT_PARAMETERS][index - 1];
}

} // namespace

namespace Kratos
{
int ConstitutiveLawUtilities::GetStateVariableIndex(const Variable<double>& rThisVariable)
{
    int index = -1;
    if (const std::string prefix{"STATE_VARIABLE_"}; rThisVariable.Name().substr(0, prefix.length()) == prefix) {
        index = std::stoi(rThisVariable.Name().substr(prefix.length()));
    }

    return index - 1;
}

void ConstitutiveLawUtilities::SetConstitutiveParameters(ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                                         Vector&       rStrainVector,
                                                         Matrix&       rConstitutiveMatrix,
                                                         const Vector& rN,
                                                         const Matrix& rGradNpT,
                                                         const Matrix& rF,
                                                         double        detF)
{
    rConstitutiveParameters.SetStrainVector(rStrainVector);
    rConstitutiveParameters.SetConstitutiveMatrix(rConstitutiveMatrix);
    rConstitutiveParameters.SetShapeFunctionsValues(rN);
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rGradNpT);
    rConstitutiveParameters.SetDeformationGradientF(rF);
    rConstitutiveParameters.SetDeterminantF(detF);
}

double ConstitutiveLawUtilities::GetCohesion(const Properties& rProperties)
{
    if (rProperties.Has(GEO_COHESION)) return rProperties[GEO_COHESION];

    try {
        return GetValueOfUMatParameter(rProperties, INDEX_OF_UMAT_C_PARAMETER);
    } catch (const std::exception& e) {
        KRATOS_ERROR
            << "ConstitutiveLawUtilities::GetCohesion failed. There is no GEO_COHESION available "
               "and attempting to get the cohesion from UMAT parameters resulted in the following "
            << e.what() << "." << std::endl;
    }
}

bool ConstitutiveLawUtilities::HasFrictionAngle(const Properties& rProperties)
{
    // Friction angle can be supplied either directly via GEO_FRICTION_ANGLE
    // or via UMAT parameters + INDEX_OF_UMAT_PHI_PARAMETER.
    return rProperties.Has(GEO_FRICTION_ANGLE) ||
           (rProperties.Has(INDEX_OF_UMAT_PHI_PARAMETER) && rProperties.Has(UMAT_PARAMETERS));
}

void ConstitutiveLawUtilities::ValidateFrictionAngle(const Properties& rProperties, IndexType ElementId)
{
    double      phi = 0.0;
    std::string phi_name;

    if (rProperties.Has(GEO_FRICTION_ANGLE)) {
        phi      = rProperties[GEO_FRICTION_ANGLE];
        phi_name = "GEO_FRICTION_ANGLE";
    } else if (rProperties.Has(INDEX_OF_UMAT_PHI_PARAMETER) && rProperties.Has(UMAT_PARAMETERS)) {
        const auto phi_index = rProperties[INDEX_OF_UMAT_PHI_PARAMETER];
        const auto number_of_umat_parameters = static_cast<int>(rProperties[UMAT_PARAMETERS].size());

        KRATOS_ERROR_IF(phi_index < 1 || phi_index > number_of_umat_parameters)
            << "Properties ( " << rProperties.Id() << ") of element ( " << ElementId
            << "): INDEX_OF_UMAT_PHI_PARAMETER (" << phi_index
            << ") is not in range [1, size of UMAT_PARAMETERS]." << std::endl;

        phi      = rProperties[UMAT_PARAMETERS][phi_index - 1];
        phi_name = "Phi";
    }

    if (phi_name.empty()) {
        KRATOS_ERROR << "Properties ( " << rProperties.Id() << ") of element ( " << ElementId
                     << ") does not have GEO_FRICTION_ANGLE nor INDEX_OF_UMAT_PHI_PARAMETER." << std::endl;
    }

    KRATOS_ERROR_IF(phi < 0.0 || phi > 90.0)
        << "Properties ( " << rProperties.Id() << ") of element ( " << ElementId << "): " << phi_name
        << " (" << phi << " degrees) should be between 0 and 90 degrees." << std::endl;
}

double ConstitutiveLawUtilities::GetFrictionAngleInDegrees(const Properties& rProperties)
{
    if (rProperties.Has(GEO_FRICTION_ANGLE)) {
        return rProperties[GEO_FRICTION_ANGLE];
    }

    try {
        return GetValueOfUMatParameter(rProperties, INDEX_OF_UMAT_PHI_PARAMETER);
    } catch (const std::exception& e) {
        KRATOS_ERROR
            << "ConstitutiveLawUtilities::GetFrictionAngleInDegrees failed. There is no "
               "GEO_FRICTION_ANGLE available and attempting to get the friction angle from UMAT "
               "parameters resulted in the following "
            << e.what() << "." << std::endl;
    }
}

double ConstitutiveLawUtilities::GetFrictionAngleInRadians(const Properties& rProperties)
{
    return MathUtils<>::DegreesToRadians(GetFrictionAngleInDegrees(rProperties));
}

Matrix ConstitutiveLawUtilities::MakeInterfaceElasticConstitutiveTensor(double      NormalStiffness,
                                                                        double      ShearStiffness,
                                                                        std::size_t TractionSize,
                                                                        std::size_t NumberOfNormalComponents)
{
    auto result = Matrix{ZeroMatrix{TractionSize, TractionSize}};
    for (auto i = std::size_t{0}; i < NumberOfNormalComponents; ++i)
        result(i, i) = NormalStiffness;
    for (auto i = NumberOfNormalComponents; i < TractionSize; ++i)
        result(i, i) = ShearStiffness;
    return result;
}

void ConstitutiveLawUtilities::CheckStrainSize(const Properties& rProperties, std::size_t ExpectedSize, std::size_t ElementId)
{
    const std::size_t strain_size = rProperties[CONSTITUTIVE_LAW]->GetStrainSize();
    KRATOS_ERROR_IF_NOT(strain_size == ExpectedSize)
        << "Wrong constitutive law is used: strain size is " << strain_size << " when it is expected to be "
        << ExpectedSize << " at element Id = " << ElementId << "." << std::endl;
}

void ConstitutiveLawUtilities::CheckHasStrainMeasure_Infinitesimal(const Properties& rProperties, std::size_t ElementId)
{
    ConstitutiveLaw::Features LawFeatures;
    rProperties[CONSTITUTIVE_LAW]->GetLawFeatures(LawFeatures);
    const auto correct_strain_measure =
        std::ranges::any_of(LawFeatures.mStrainMeasures, [](auto& strain_measure) {
        return strain_measure == ConstitutiveLaw::StrainMeasure_Infinitesimal;
    });

    KRATOS_ERROR_IF_NOT(correct_strain_measure)
        << "Constitutive law is not compatible with the strain type "
           "StrainMeasure_Infinitesimal at element "
        << ElementId << "." << std::endl;
}

double ConstitutiveLawUtilities::CalculateK0NCFromFrictionAngleInRadians(double FrictionAngleInRadians)
{
    return 1.0 - std::sin(FrictionAngleInRadians);
}

double ConstitutiveLawUtilities::CalculateUndrainedYoungsModulus(const Properties& rProperties, double UndrainedPoissonsRatio)
{
    const auto denominator = 1.0 + rProperties[POISSON_RATIO];
    return rProperties[YOUNG_MODULUS] * (1.0 + UndrainedPoissonsRatio) / denominator;
}

double ConstitutiveLawUtilities::GetOrCalculateUndrainedPoissonsRatio(const Properties& rProperties)
{
    const auto result = rProperties.Has(GEO_POISSON_UNDRAINED)
                            ? rProperties[GEO_POISSON_UNDRAINED]
                            : CalculateUndrainedPoissonsRatio(rProperties);

    if (constexpr auto max_poisson_ratio = 0.495; result > max_poisson_ratio) {
        KRATOS_WARNING("Clamping undrained Poisson ratio from ")
            << result << " to " << max_poisson_ratio << "." << std::endl;
        return max_poisson_ratio;
    }

    return result;
}

double ConstitutiveLawUtilities::CalculateUndrainedPoissonsRatio(const Properties& rProperties)
{
    const auto skempton_b       = GetOrCalculateSkemptonB(rProperties);
    const auto biot_coefficient = rProperties[BIOT_COEFFICIENT];
    const auto poissons_ratio   = rProperties[POISSON_RATIO];
    const auto denominator = 3.0 - biot_coefficient * skempton_b * (1.0 - 2.0 * poissons_ratio);
    KRATOS_ERROR_IF(denominator <= std::numeric_limits<double>::epsilon())
        << "Non-physical values: denominator < epsilon." << std::endl;
    return (3.0 * poissons_ratio + biot_coefficient * skempton_b * (1.0 - 2.0 * poissons_ratio)) / denominator;
}

double ConstitutiveLawUtilities::GetOrCalculateSkemptonB(const Properties& rProperties)
{
    if (rProperties.Has(GEO_SKEMPTON_B)) {
        return rProperties[GEO_SKEMPTON_B];
    }

    const auto k_f = rProperties[BULK_MODULUS_FLUID];
    const auto k_s = rProperties[BULK_MODULUS_SOLID]; // assuming the bulk stiffness of the skeleton equals the solid bulk stiffness
    const auto porosity         = rProperties[POROSITY];
    const auto biot_coefficient = rProperties[BIOT_COEFFICIENT];
    const auto denominator = biot_coefficient + porosity * ((k_s / k_f) + biot_coefficient - 1.0);
    KRATOS_ERROR_IF(denominator <= std::numeric_limits<double>::epsilon())
        << "Non-physical values: denominator < epsilon." << denominator << std::endl;

    const auto result = biot_coefficient / denominator;
    KRATOS_ERROR_IF(result < -std::numeric_limits<double>::epsilon() || result > 1.0 + std::numeric_limits<double>::epsilon())
        << "Calculated Skempton B (" << result << ") is out of range [0,1]." << std::endl;

    return std::clamp(result, 0.0, 1.0);
}

std::pair<double, double> ConstitutiveLawUtilities::GetOrCalculateElasticProperties(const Properties& rProperties,
                                                                                    bool Undrained)
{
    if (Undrained) {
        const auto nu = GetOrCalculateUndrainedPoissonsRatio(rProperties);
        const auto E  = CalculateUndrainedYoungsModulus(rProperties, nu);
        return {E, nu};
    }
    return {rProperties[YOUNG_MODULUS], rProperties[POISSON_RATIO]};
}

Matrix ConstitutiveLawUtilities::MakeContinuumElasticConstitutiveTensor(double      YoungsModulus,
                                                                        double      PoissonsRatio,
                                                                        std::size_t StrainSize,
                                                                        std::size_t NumberOfNormalComponents)
{
    const auto denominator = (1.0 + PoissonsRatio) * (1.0 - 2.0 * PoissonsRatio);
    KRATOS_ERROR_IF(denominator <= std::numeric_limits<double>::epsilon())
        << "PoissonsRatio of " << PoissonsRatio << " leads to a nearly zero denominator" << std::endl;
    const auto prefactor          = YoungsModulus / denominator;
    const auto diagonal_value     = (1.0 - PoissonsRatio) * prefactor;
    const auto off_diagonal_value = PoissonsRatio * prefactor;

    auto result = Matrix{ZeroMatrix{StrainSize, StrainSize}};
    for (auto i = std::size_t{0}; i < NumberOfNormalComponents; ++i) {
        for (auto j = std::size_t{0}; j < NumberOfNormalComponents; ++j) {
            result(i, j) = i == j ? diagonal_value : off_diagonal_value;
        }
    }
    const auto shear_modulus = YoungsModulus / (2.0 * (1.0 + PoissonsRatio));
    for (auto i = NumberOfNormalComponents; i < StrainSize; ++i) {
        result(i, i) = shear_modulus;
    }
    return result;
}

DrainageType ConstitutiveLawUtilities::StringToDrainageType(const std::string& rDrainageTypeName)
{
    using enum DrainageType;
    static const std::map<std::string, DrainageType, std::less<>> drainage_type_map = {
        {"drained"s, DRAINED},
        {"fully_coupled"s, FULLY_COUPLED},
        {"undrained"s, UNDRAINED},
        {"constant_pw_field"s, CONSTANT_WATER_PRESSURE}};
    return drainage_type_map.at(GeoStringUtilities::ToLower(rDrainageTypeName));
}

bool ConstitutiveLawUtilities::IsConstantWaterPressure(const Properties& rProperties)
{
    return ConstitutiveLawUtilities::StringToDrainageType(rProperties[GEO_DRAINAGE_TYPE]) ==
           DrainageType::CONSTANT_WATER_PRESSURE;
}

void ConstitutiveLawUtilities::ReplaceIgnoreUndrainedByDrainageType(Properties& rProperties)
{
    if (!rProperties.Has(IGNORE_UNDRAINED)) return;

    KRATOS_WARNING("DEPRECATION")
        << "Use of IGNORE_UNDRAINED is deprecated, please change your input to "
           "GEO_DRAINAGE_TYPE"
        << std::endl;
    if (!rProperties.Has(GEO_DRAINAGE_TYPE))
        rProperties[GEO_DRAINAGE_TYPE] = rProperties[IGNORE_UNDRAINED] ? "constant_pw_field"s : "fully_coupled"s;
    rProperties.Erase(IGNORE_UNDRAINED);
}

double ConstitutiveLawUtilities::CalculateExcessPorePressureIncrement(const Properties& rProperties,
                                                                      double VolumetricStrainIncrement)
{
    const auto bulk_modulus_fluid = rProperties[BULK_MODULUS_FLUID];
    const auto bulk_modulus_solid = rProperties[BULK_MODULUS_SOLID];
    const auto porosity           = rProperties[POROSITY];
    const auto biot_coefficient   = rProperties[BIOT_COEFFICIENT];

    const double denominator = porosity / bulk_modulus_fluid + (biot_coefficient - porosity) / bulk_modulus_solid;
    KRATOS_ERROR_IF(std::abs(denominator) <= std::numeric_limits<double>::epsilon())
        << "Non-physical values: denominator < epsilon for property Id of " << rProperties.Id()
        << "." << std::endl;
    return biot_coefficient * VolumetricStrainIncrement / denominator;
}

} // namespace Kratos