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
//  Main authors:    Mohamed Nabi,
//                   Wijtze Pieter Kikstra
//

// Application includes
#include "custom_constitutive/mohr_coulomb_with_tension_cutoff.h"
#include "custom_constitutive/constitutive_law_dimension.h"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

#include <cmath>
#include <type_traits>

namespace
{

using namespace Kratos;

void CheckProperty(const Properties&       rMaterialProperties,
                   const Variable<double>& rVariable,
                   std::optional<double>   MaxValue = std::nullopt)
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(rVariable))
        << rVariable.Name() << " is not defined for property " << rMaterialProperties.Id() << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[rVariable] < 0.0 ||
                    (MaxValue.has_value() && rMaterialProperties[rVariable] > MaxValue.value()))
        << "value of " << rVariable.Name() << " for property " << rMaterialProperties.Id()
        << " is out of range: " << rMaterialProperties[rVariable] << " is not in [0.0, "
        << (MaxValue ? std::to_string(*MaxValue) + "]" : "->") << std::endl;
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

bool IsStressAtTensionApexReturnZone(const Vector& rTrialSigmaTau, double TensileStrength, double Apex)
{
    return TensileStrength < Apex && rTrialSigmaTau[0] - rTrialSigmaTau[1] - TensileStrength > 0.0;
}

bool IsStressAtTensionCutoffReturnZone(const Vector& rTrialSigmaTau, double TensileStrength, double Apex, const Vector& rCornerPoint)
{
    return TensileStrength < Apex &&
           rCornerPoint[1] - rTrialSigmaTau[1] - rCornerPoint[0] + rTrialSigmaTau[0] > 0.0;
}

bool IsStressAtCornerReturnZone(const Vector& rTrialSigmaTau, double DilatancyAngle, const Vector& rCornerPoint)
{
    return rTrialSigmaTau[0] - rCornerPoint[0] - (rTrialSigmaTau[1] - rCornerPoint[1]) * std::sin(DilatancyAngle) >= 0.0;
}

} // namespace

namespace Kratos
{
MohrCoulombWithTensionCutOff::MohrCoulombWithTensionCutOff(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mStressVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStressVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStrainVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize()))
{
}

ConstitutiveLaw::Pointer MohrCoulombWithTensionCutOff::Clone() const
{
    auto p_result = std::make_shared<MohrCoulombWithTensionCutOff>(mpConstitutiveDimension->Clone());
    p_result->mStressVector          = mStressVector;
    p_result->mStressVectorFinalized = mStressVectorFinalized;
    p_result->mStrainVectorFinalized = mStrainVectorFinalized;
    p_result->mCoulombYieldSurface   = mCoulombYieldSurface;
    p_result->mTensionCutOff         = mTensionCutOff;
    return p_result;
}

Vector& MohrCoulombWithTensionCutOff::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == CAUCHY_STRESS_VECTOR) {
        rValue = mStressVector;
    } else {
        rValue = ConstitutiveLaw::GetValue(rThisVariable, rValue);
    }
    return rValue;
}

void MohrCoulombWithTensionCutOff::SetValue(const Variable<Vector>& rVariable,
                                            const Vector&           rValue,
                                            const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == CAUCHY_STRESS_VECTOR) {
        mStressVector = rValue;
    } else {
        KRATOS_ERROR << "Can't set value of " << rVariable.Name() << ": unsupported variable\n";
    }
}

SizeType MohrCoulombWithTensionCutOff::WorkingSpaceDimension()
{
    return mpConstitutiveDimension->GetDimension();
}

int MohrCoulombWithTensionCutOff::Check(const Properties&   rMaterialProperties,
                                        const GeometryType& rElementGeometry,
                                        const ProcessInfo&  rCurrentProcessInfo) const
{
    const auto result = ConstitutiveLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    CheckProperty(rMaterialProperties, GEO_COHESION);
    CheckProperty(rMaterialProperties, GEO_FRICTION_ANGLE);
    CheckProperty(rMaterialProperties, GEO_DILATANCY_ANGLE, rMaterialProperties[GEO_FRICTION_ANGLE]);
    CheckProperty(rMaterialProperties, GEO_TENSILE_STRENGTH,
                  rMaterialProperties[GEO_COHESION] /
                      std::tan(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE])));
    CheckProperty(rMaterialProperties, YOUNG_MODULUS);
    CheckProperty(rMaterialProperties, POISSON_RATIO, 0.5);
    return result;
}

ConstitutiveLaw::StressMeasure MohrCoulombWithTensionCutOff::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

SizeType MohrCoulombWithTensionCutOff::GetStrainSize() const
{
    return mpConstitutiveDimension->GetStrainSize();
}

ConstitutiveLaw::StrainMeasure MohrCoulombWithTensionCutOff::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

bool MohrCoulombWithTensionCutOff::IsIncremental() { return true; }

bool MohrCoulombWithTensionCutOff::RequiresInitializeMaterialResponse() { return true; }

void MohrCoulombWithTensionCutOff::InitializeMaterial(const Properties& rMaterialProperties,
                                                      const Geometry<Node>&,
                                                      const Vector&)
{
    mCoulombYieldSurface =
        CoulombYieldSurface(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE]),
                            rMaterialProperties[GEO_COHESION],
                            MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_DILATANCY_ANGLE]));
    mTensionCutOff = TensionCutoff(rMaterialProperties[GEO_TENSILE_STRENGTH]);
}

void MohrCoulombWithTensionCutOff::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    if (!mIsModelInitialized) {
        mStressVectorFinalized = rValues.GetStressVector();
        mStrainVectorFinalized = rValues.GetStrainVector();
        mIsModelInitialized    = true;
    }
}

void MohrCoulombWithTensionCutOff::GetLawFeatures(Features& rFeatures)
{
    auto options = Flags{};
    options.Set(mpConstitutiveDimension->GetSpatialType());
    options.Set(ConstitutiveLaw::INFINITESIMAL_STRAINS);
    options.Set(ConstitutiveLaw::ISOTROPIC);
    rFeatures.SetOptions(options);

    rFeatures.SetStrainMeasure(StrainMeasure_Infinitesimal);
    rFeatures.SetStrainSize(GetStrainSize());
    rFeatures.SetSpaceDimension(WorkingSpaceDimension());
}

void MohrCoulombWithTensionCutOff::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters)
{
    const auto& r_prop = rParameters.GetMaterialProperties();

    if (rParameters.GetOptions().Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        rParameters.GetConstitutiveMatrix() =
            mpConstitutiveDimension->CalculateElasticMatrix(r_prop[YOUNG_MODULUS], r_prop[POISSON_RATIO]);
    }

    const auto trail_stress_vector = CalculateTrialStressVector(
        rParameters.GetStrainVector(), r_prop[YOUNG_MODULUS], r_prop[POISSON_RATIO]);

    Vector principal_trial_stress_vector;
    Matrix rotation_matrix;
    StressStrainUtilities::CalculatePrincipalStresses(
        trail_stress_vector, principal_trial_stress_vector, rotation_matrix);

    while (!IsAdmissiblePrincipalStressState(principal_trial_stress_vector)) {
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
                mTensionCutOff.DerivativeOfFlowFunction(principal_trial_stress_vector),
                r_prop[GEO_TENSILE_STRENGTH]);
        } else if (IsStressAtCornerReturnZone(
                       trial_sigma_tau, MathUtils<>::DegreesToRadians(r_prop[GEO_DILATANCY_ANGLE]), corner_point)) {
            principal_trial_stress_vector =
                ReturnStressAtCornerReturnZone(principal_trial_stress_vector, corner_point);
        } else {
            // Regular failure region
            principal_trial_stress_vector = ReturnStressAtRegularFailureZone(
                principal_trial_stress_vector,
                mCoulombYieldSurface.DerivativeOfFlowFunction(principal_trial_stress_vector),
                MathUtils<>::DegreesToRadians(r_prop[GEO_FRICTION_ANGLE]), r_prop[GEO_COHESION]);
        }

        StressStrainUtilities::ReorderEigenValuesAndVectors(principal_trial_stress_vector, rotation_matrix);
    }

    mStressVector = StressStrainUtilities::RotatePrincipalStresses(
        principal_trial_stress_vector, rotation_matrix, mpConstitutiveDimension->GetStrainSize());

    rParameters.GetStressVector() = mStressVector;
}

bool MohrCoulombWithTensionCutOff::IsAdmissiblePrincipalStressState(const Vector& rPrincipalStresses) const
{
    constexpr auto tolerance          = 1.0e-10;
    const auto coulomb_yield_function = mCoulombYieldSurface.YieldFunctionValue(rPrincipalStresses);
    const auto tension_yield_function = mTensionCutOff.YieldFunctionValue(rPrincipalStresses);
    const auto coulomb_tolerance      = tolerance * (1.0 + std::abs(coulomb_yield_function));
    const auto tension_tolerance      = tolerance * (1.0 + std::abs(tension_yield_function));
    return coulomb_yield_function < coulomb_tolerance && tension_yield_function < tension_tolerance;
}

Vector MohrCoulombWithTensionCutOff::CalculateTrialStressVector(const Vector& rStrainVector,
                                                                double        YoungsModulus,
                                                                double        PoissonsRatio) const
{
    return mStressVectorFinalized +
           prod(mpConstitutiveDimension->CalculateElasticMatrix(YoungsModulus, PoissonsRatio),
                rStrainVector - mStrainVectorFinalized);
}

void MohrCoulombWithTensionCutOff::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mStrainVectorFinalized = rValues.GetStrainVector();
    mStressVectorFinalized = mStressVector;
}

void MohrCoulombWithTensionCutOff::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.save("ConstitutiveLawDimension", mpConstitutiveDimension);
    rSerializer.save("StressVector", mStressVector);
    rSerializer.save("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.save("CoulombYieldSurface", mCoulombYieldSurface);
    rSerializer.save("TensionCutOff", mTensionCutOff);
    rSerializer.save("IsModelInitialized", mIsModelInitialized);
}

void MohrCoulombWithTensionCutOff::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.load("ConstitutiveLawDimension", mpConstitutiveDimension);
    rSerializer.load("StressVector", mStressVector);
    rSerializer.load("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
    rSerializer.load("CoulombYieldSurface", mCoulombYieldSurface);
    rSerializer.load("TensionCutOff", mTensionCutOff);
    rSerializer.load("IsModelInitialized", mIsModelInitialized);
}

static_assert(!std::is_copy_constructible_v<MohrCoulombWithTensionCutOff>);
static_assert(!std::is_copy_assignable_v<MohrCoulombWithTensionCutOff>);
static_assert(std::is_move_constructible_v<MohrCoulombWithTensionCutOff>);
static_assert(std::is_move_assignable_v<MohrCoulombWithTensionCutOff>);

} // Namespace Kratos