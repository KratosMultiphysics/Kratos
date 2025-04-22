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
#include "geo_mechanics_application_variables.h"

namespace
{

using namespace Kratos;

double GetValueOfUMatParameter(const Properties& rProperties, const Variable<int>& rIndexVariable)
{
    KRATOS_ERROR_IF_NOT(rProperties.Has(UMAT_PARAMETERS))
        << "Material " << rProperties.Id() << " does not have UMAT_PARAMETERS\n";

    KRATOS_ERROR_IF_NOT(rProperties.Has(rIndexVariable))
        << "Material " << rProperties.Id() << " does not have " << rIndexVariable.Name() << "\n";

    const auto index = rProperties[rIndexVariable]; // 1-based index
    KRATOS_DEBUG_ERROR_IF(index < 1 ||
                          static_cast<std::size_t>(index) > rProperties[UMAT_PARAMETERS].size())
        << "Got out-of-bounds " << rIndexVariable.Name() << " (material ID: " << rProperties.Id()
        << "): " << index << " is not in range [1, " << rProperties[UMAT_PARAMETERS].size() << "]\n";

    return rProperties[UMAT_PARAMETERS][index - 1];
}



Vector TransformPrincipalStressesToSigmaAndTau(const Vector& rPrincipalStresses)
{
    auto result = Vector(2);
    result[0]   = 0.5 * (rPrincipalStresses[0] + rPrincipalStresses[2]);
    result[1]   = 0.5 * (rPrincipalStresses[0] - rPrincipalStresses[2]);
    return result;
}

double CalculateApex(double FrictionAngle, double Cohesion)
{
    return Cohesion / std::tan(FrictionAngle);
}

Vector CalculateCornerPoint(double FrictionAngle, double Cohesion, double TensileStrength)
{
    // Check whether the tension cut-off lies beyond the apex
    auto result = Vector{ZeroVector(2)};
    result[0]   = CalculateApex(FrictionAngle, Cohesion);
    if (TensileStrength > result[0]) return result;

    result[0] = (TensileStrength - Cohesion * std::cos(FrictionAngle)) / (1.0 - std::sin(FrictionAngle));
    result[1] = (Cohesion * std::cos(FrictionAngle) - TensileStrength * std::sin(FrictionAngle)) /
                (1.0 - std::sin(FrictionAngle));
    return result;
}

Vector ReturnStressAtTensionApexReturnZone(const Vector& rPrincipalTrialStressVector, double TensileStrength)
{
    auto result = rPrincipalTrialStressVector;
    result[0]   = TensileStrength;
    result[2]   = result[0];
    return result;
}

Vector ReturnStressAtTensionCutoffReturnZone(const Vector& rPrincipalTrialStressVector,
                                             const Vector& rDerivativeOfFlowFunction,
                                             double        TensileStrength)
{
    const auto lambda = (TensileStrength - rPrincipalTrialStressVector[0]) / rDerivativeOfFlowFunction[0];
    return rPrincipalTrialStressVector + lambda * rDerivativeOfFlowFunction;
}

Vector ReturnStressAtCornerReturnZone(const Vector& rPrincipalTrialStressVector, const Vector& rCornerPoint)
{
    auto result = rPrincipalTrialStressVector;
    result[0]   = rCornerPoint[0] + rCornerPoint[1];
    result[2]   = rCornerPoint[0] - rCornerPoint[1];
    return result;
}

Vector ReturnStressAtRegularFailureZone(const Vector& rPrincipalTrialStressVector,
                                        const Vector& rDerivativeOfFlowFunction,
                                        double        FrictionAngle,
                                        double        Cohesion)
{
    const auto cof1 = (1.0 + std::sin(FrictionAngle)) / (1.0 - std::sin(FrictionAngle));
    const auto cof2 = 2.0 * Cohesion * std::cos(FrictionAngle) / (1.0 - std::sin(FrictionAngle));
    const auto numerator = cof1 * rDerivativeOfFlowFunction[0] - rDerivativeOfFlowFunction[2];
    const auto lambda =
        (rPrincipalTrialStressVector[2] + cof2 - rPrincipalTrialStressVector[0] * cof1) / numerator;
    return rPrincipalTrialStressVector + lambda * rDerivativeOfFlowFunction;
}

bool IsStressAtTensionApexReturnZone(const Vector& rTrialSigmaTau,
                                                                   double        TensileStrength,
                                                                   double        Apex)
{
    return TensileStrength < Apex && rTrialSigmaTau[0] - rTrialSigmaTau[1] - TensileStrength > 0.0;
}

bool IsStressAtTensionCutoffReturnZone(const Vector& rTrialSigmaTau,
                                                                     double        TensileStrength,
                                                                     double        Apex,
                                                                     const Vector& rCornerPoint)
{
    return TensileStrength < Apex &&
           rCornerPoint[1] - rTrialSigmaTau[1] - rCornerPoint[0] + rTrialSigmaTau[0] > 0.0;
}

bool IsStressAtCornerReturnZone(const Vector& rTrialSigmaTau,
                                                              double        DilatancyAngle,
                                                              const Vector& rCornerPoint)
{
    return rTrialSigmaTau[0] - rCornerPoint[0] - (rTrialSigmaTau[1] - rCornerPoint[1]) * std::sin(DilatancyAngle) >= 0.0;
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
    return rProperties.Has(GEO_COHESION) ? rProperties[GEO_COHESION]
                                         : GetValueOfUMatParameter(rProperties, INDEX_OF_UMAT_C_PARAMETER);
}

double ConstitutiveLawUtilities::GetFrictionAngleInDegrees(const Properties& rProperties)
{
    return rProperties.Has(GEO_FRICTION_ANGLE)
               ? rProperties[GEO_FRICTION_ANGLE]
               : GetValueOfUMatParameter(rProperties, INDEX_OF_UMAT_PHI_PARAMETER);
}


Vector ConstitutiveLawUtilities::MapStressesInMorhCoulomb(const Properties& r_prop,
                                                          Vector& principal_trial_stress_vector,
                                                          const CoulombYieldSurface& rCoulombYieldSurface,
                                                          const TensionCutoff& rTensionCutOff)
{
        const auto apex = CalculateApex(MathUtils<>::DegreesToRadians(r_prop[GEO_FRICTION_ANGLE]),
                                        r_prop[GEO_COHESION]);
        const auto corner_point =
            CalculateCornerPoint(MathUtils<>::DegreesToRadians(r_prop[GEO_FRICTION_ANGLE]),
                                 r_prop[GEO_COHESION], r_prop[GEO_TENSILE_STRENGTH]);

        if (const auto trial_sigma_tau = TransformPrincipalStressesToSigmaAndTau(principal_trial_stress_vector);
            IsStressAtTensionApexReturnZone(trial_sigma_tau, r_prop[GEO_TENSILE_STRENGTH], apex)) {
            principal_trial_stress_vector = ReturnStressAtTensionApexReturnZone(
                principal_trial_stress_vector, r_prop[GEO_TENSILE_STRENGTH]);
        } else if (IsStressAtTensionCutoffReturnZone(trial_sigma_tau, r_prop[GEO_TENSILE_STRENGTH],
                                                     apex, corner_point)) {
            principal_trial_stress_vector = ReturnStressAtTensionCutoffReturnZone(
                principal_trial_stress_vector,
                rTensionCutOff.DerivativeOfFlowFunction(principal_trial_stress_vector),
                r_prop[GEO_TENSILE_STRENGTH]);
        } else if (IsStressAtCornerReturnZone(
                       trial_sigma_tau, MathUtils<>::DegreesToRadians(r_prop[GEO_DILATANCY_ANGLE]), corner_point)) {
            principal_trial_stress_vector =
                ReturnStressAtCornerReturnZone(principal_trial_stress_vector, corner_point);
        } else {
            // Regular failure region
            principal_trial_stress_vector = ReturnStressAtRegularFailureZone(
                principal_trial_stress_vector,
                rCoulombYieldSurface.DerivativeOfFlowFunction(principal_trial_stress_vector),
                MathUtils<>::DegreesToRadians(r_prop[GEO_FRICTION_ANGLE]), r_prop[GEO_COHESION]);
        }
    return principal_trial_stress_vector;
}


bool ConstitutiveLawUtilities::IsAdmissiblePrincipalStressState(const Vector& rPrincipalStresses,
                                                                const CoulombYieldSurface& rCoulombYieldSurface,
                                                                const TensionCutoff& rTensionCutOff)
{
    constexpr auto tolerance          = 1.0e-10;
    const auto coulomb_yield_function = rCoulombYieldSurface.YieldFunctionValue(rPrincipalStresses);
    const auto tension_yield_function = rTensionCutOff.YieldFunctionValue(rPrincipalStresses);
    const auto coulomb_tolerance      = tolerance * (1.0 + std::abs(coulomb_yield_function));
    const auto tension_tolerance      = tolerance * (1.0 + std::abs(tension_yield_function));
    return coulomb_yield_function < coulomb_tolerance && tension_yield_function < tension_tolerance;
}

} // namespace Kratos
