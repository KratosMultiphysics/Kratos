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
#include "includes/serializer.h"

namespace
{

using namespace Kratos;

Vector TransformPrincipalStressesToSigmaAndTau(const Vector& rPrincipalStresses)
{
    auto result = Vector(2);
    result[0]   = 0.5 * (rPrincipalStresses(0) + rPrincipalStresses(2));
    result[1]   = 0.5 * (rPrincipalStresses(0) - rPrincipalStresses(2));
    return result;
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
    p_result->mCoulombYieldSurface = mCoulombYieldSurface;
    p_result->mTensionCutOff = mTensionCutOff;
    return p_result;
}

Vector& MohrCoulombWithTensionCutOff::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == CAUCHY_STRESS_VECTOR) {
        rValue = mStressVector;
    } else {
        KRATOS_ERROR << "Can't get value of " << rThisVariable.Name() << ": unsupported variable\n";
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
    CheckProperty(rMaterialProperties, GEO_DILATANCY_ANGLE);
    CheckProperty(rMaterialProperties, GEO_TENSILE_STRENGTH);
    CheckProperty(rMaterialProperties, YOUNG_MODULUS);
    CheckProperty(rMaterialProperties, POISSON_RATIO);
    return result;
}

void MohrCoulombWithTensionCutOff::CheckProperty(const Properties& rMaterialProperties,
                                                 const Kratos::Variable<double>& rVariable) const
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(rVariable) || rMaterialProperties[rVariable] < 0.0)
        << rVariable.Name()
        << " is not defined or has an invalid value for property: " << rMaterialProperties.Id()
        << std::endl;
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
    const auto friction_angle_in_rad = MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE]);
    const auto cohesion = rMaterialProperties[GEO_COHESION];
    const auto dilatancy_angle_in_rad = MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_DILATANCY_ANGLE]);
    const auto tensile_strength = rMaterialProperties[GEO_TENSILE_STRENGTH];

    mCoulombYieldSurface = CoulombYieldSurface(friction_angle_in_rad, cohesion, dilatancy_angle_in_rad);
    mTensionCutOff       = TensionCutoff(tensile_strength);
}

void MohrCoulombWithTensionCutOff::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters)
{
    const auto& r_prop = rParameters.GetMaterialProperties();

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

        if (IsStressAtAxialZone(principal_trial_stress_vector, r_prop[GEO_TENSILE_STRENGTH], apex, corner_point)) {
            principal_trial_stress_vector =
                ReturnStressAtAxialZone(principal_trial_stress_vector, r_prop[GEO_TENSILE_STRENGTH]);
        } else if (IsStressAtCornerReturnZone(principal_trial_stress_vector,
                                              MathUtils<>::DegreesToRadians(r_prop[GEO_DILATANCY_ANGLE]),
                                              corner_point)) {
            principal_trial_stress_vector =
                ReturnStressAtCornerReturnZone(principal_trial_stress_vector, corner_point);
        } else {
            // Regular failure region
            principal_trial_stress_vector = ReturnStressAtRegularFailureZone(
                principal_trial_stress_vector,
                MathUtils<>::DegreesToRadians(r_prop[GEO_FRICTION_ANGLE]), r_prop[GEO_COHESION]);
        }

        StressStrainUtilities::ReorderEigenValuesAndVectors(principal_trial_stress_vector, rotation_matrix);
    }

    mStressVector = StressStrainUtilities::RotatePrincipalStresses(
        principal_trial_stress_vector, rotation_matrix, mpConstitutiveDimension->GetStrainSize());
}

Vector MohrCoulombWithTensionCutOff::ReturnStressAtAxialZone(const Vector& rPrincipalTrialStressVector,
                                                             double TensileStrength) const
{
    Vector result = rPrincipalTrialStressVector;
    result[0]     = TensileStrength;
    result[2] = TensileStrength - rPrincipalTrialStressVector(0) + rPrincipalTrialStressVector(2);
    return result;
}

Vector MohrCoulombWithTensionCutOff::ReturnStressAtCornerReturnZone(const Vector& rPrincipalTrialStressVector,
                                                                    const Vector& rCornerPoint) const
{
    Vector result = rPrincipalTrialStressVector;
    result[0]     = rCornerPoint[0] + rCornerPoint[1];
    result[2]     = rCornerPoint[0] - rCornerPoint[1];
    return result;
}

Vector MohrCoulombWithTensionCutOff::ReturnStressAtRegularFailureZone(const Vector& rPrincipalTrialStressVector,
                                                                      double FrictionAngle,
                                                                      double Cohesion) const
{
    const auto flow_derivative = mCoulombYieldSurface.DerivativeOfFlowFunction(rPrincipalTrialStressVector);
    const auto cof1 = (1.0 + std::sin(FrictionAngle)) / (1.0 - std::sin(FrictionAngle));
    const auto cof2 = 2.0 * Cohesion * std::cos(FrictionAngle) / (1.0 - std::sin(FrictionAngle));
    const auto numerator  = cof1 * flow_derivative[0] - flow_derivative[2];
    const auto lambda = (rPrincipalTrialStressVector[2] + cof2 - rPrincipalTrialStressVector[0] * cof1) / numerator;
    return rPrincipalTrialStressVector + lambda * flow_derivative;;
}

bool MohrCoulombWithTensionCutOff::IsAdmissiblePrincipalStressState(const Vector& rPrincipalStresses) const
{
    return mCoulombYieldSurface.YieldFunctionValue(rPrincipalStresses) <= 0.0 &&
           mTensionCutOff.YieldFunctionValue(rPrincipalStresses) <= 0.0;
}

bool MohrCoulombWithTensionCutOff::IsStressAtAxialZone(const Vector& rPrincipalTrialStresses,
                                                       double        TensileStrength,
                                                       double        Apex,
                                                       const Vector& rCornerPoint) const
{
    const auto trial_tau = TransformPrincipalStressesToSigmaAndTau(rPrincipalTrialStresses)[1];
    return TensileStrength < Apex && trial_tau <= rCornerPoint[1] &&
           mTensionCutOff.YieldFunctionValue(rPrincipalTrialStresses) >= 0.0;
}

bool MohrCoulombWithTensionCutOff::IsStressAtCornerReturnZone(const Vector& rPrincipalTrialStresses,
                                                              double        DilatancyAngle,
                                                              const Vector& rCornerPoint) const
{
    const auto trial_sigma_tau = TransformPrincipalStressesToSigmaAndTau(rPrincipalTrialStresses);
    return trial_sigma_tau[0] - rCornerPoint[0] - (trial_sigma_tau[1] - rCornerPoint[1]) * std::sin(DilatancyAngle) >= 0.0;
}

double MohrCoulombWithTensionCutOff::CalculateApex(double FrictionAngle, double Cohesion) const
{
    return Cohesion / std::tan(FrictionAngle);
}

Vector MohrCoulombWithTensionCutOff::CalculateCornerPoint(double FrictionAngle, double Cohesion, double TensileStrength) const
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
    rSerializer.save("StressVector", mStressVector);
    rSerializer.save("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.save("StrainVectorFinalized", mStrainVectorFinalized);
    //rSerializer.save("CoulombYieldSurface", mCoulombYieldSurface);
    //rSerializer.save("TensionCutOff", mTensionCutOff);
}

void MohrCoulombWithTensionCutOff::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.load("StressVector", mStressVector);
    rSerializer.load("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
    //rSerializer.load("CoulombYieldSurface", mCoulombYieldSurface);
    //rSerializer.load("TensionCutOff", mTensionCutOff);
}
} // Namespace Kratos