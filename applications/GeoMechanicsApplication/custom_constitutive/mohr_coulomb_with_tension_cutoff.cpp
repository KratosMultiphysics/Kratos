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
#include "custom_constitutive/mohr_coulomb_with_tension_cutoff.hpp"
#include "custom_constitutive/constitutive_law_dimension.h"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"

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

void MohrCoulombWithTensionCutOff::SetValue(const Variable<Vector >& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
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
    ConstitutiveLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    CheckProperty(rMaterialProperties, GEO_COHESION);
    CheckProperty(rMaterialProperties, GEO_FRICTION_ANGLE);
    CheckProperty(rMaterialProperties, GEO_DILATANCY_ANGLE);
    CheckProperty(rMaterialProperties, GEO_TENSILE_STRENGTH);
    CheckProperty(rMaterialProperties, YOUNG_MODULUS);
    CheckProperty(rMaterialProperties, POISSON_RATIO);
    return 0;
}

void MohrCoulombWithTensionCutOff::CheckProperty(const Properties& rMaterialProperties,
                                                 const Kratos::Variable<double>& rVariable) const
{
    constexpr auto min_value = 0.0;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(rVariable) || rMaterialProperties[rVariable] < min_value)
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

void MohrCoulombWithTensionCutOff::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters)
{
    const Properties& r_prop         = rParameters.GetMaterialProperties();
    const auto        friction_angle_in_rad = MathUtils<>::DegreesToRadians(r_prop[GEO_FRICTION_ANGLE]);
    const auto        cohesion       = r_prop[GEO_COHESION];
    const auto        dilatancy_angle = MathUtils<>::DegreesToRadians(r_prop[GEO_DILATANCY_ANGLE]);
    const auto        tension_cutoff = r_prop[GEO_TENSILE_STRENGTH];

    mCoulombYieldSurface = CoulombYieldSurface(friction_angle_in_rad, cohesion, dilatancy_angle);
    mTensionCutOff       = TensionCutoff(tension_cutoff);

    Vector trail_stress_vector;
    Vector strain_vector = rParameters.GetStrainVector();

    this->CalculateTrialStressVector(strain_vector, trail_stress_vector, rParameters);
    Vector principal_trial_stress_vector;
    Matrix eigenvectors_matrix;
    StressStrainUtilities::CalculatePrincipalStresses(
        trail_stress_vector, principal_trial_stress_vector, eigenvectors_matrix);
    Matrix rotation_matrix = StressStrainUtilities::CalculateRotationMatrix(eigenvectors_matrix);

    double coulomb = mCoulombYieldSurface.YieldFunctionValue(principal_trial_stress_vector);
    double cutoff  = mTensionCutOff.YieldFunctionValue(principal_trial_stress_vector);

    // Elastic region
    if (coulomb <= 0.0 && cutoff <= 0.0) {
        mStressVector = trail_stress_vector;
        return;
    }

    // Tension cut-off return
    double trail_sigma = 0.5 * (principal_trial_stress_vector(0) + principal_trial_stress_vector(2));
    double trial_tau = 0.5 * (principal_trial_stress_vector(0) - principal_trial_stress_vector(2));
    double apex      = this->CalculateApex(friction_angle_in_rad, cohesion);
    Vector corner_point = this->CalculateCornerPoint(friction_angle_in_rad, cohesion, tension_cutoff);

    if (tension_cutoff < apex && trial_tau <= corner_point(1) && cutoff >= 0.0) {
        Vector modified_principal = this->ReturnStressAtAxialZone(principal_trial_stress_vector, tension_cutoff);
        StressStrainUtilities::ReorderEigenValuesAndVectors(modified_principal, rotation_matrix);
        mStressVector = this->RotatePrincipalStresses(modified_principal, rotation_matrix);
        return;
    }

    // Zone of tensile corner return
    double derivative_flow_function =
        (trial_tau - corner_point(1)) * std::sin(dilatancy_angle) - (trail_sigma - corner_point(0));
    if (derivative_flow_function <= 0.0) {
        Vector modified_principal =
            this->ReturnStressAtCornerReturnZone(principal_trial_stress_vector, corner_point);
        StressStrainUtilities::ReorderEigenValuesAndVectors(modified_principal, rotation_matrix);
        mStressVector = this->RotatePrincipalStresses(modified_principal, rotation_matrix);
        return;
    }

    // Regular failure region
    Vector modified_principal = this->ReturnStressAtRegularFailureZone(
        principal_trial_stress_vector, mCoulombYieldSurface, friction_angle_in_rad, cohesion);
    StressStrainUtilities::ReorderEigenValuesAndVectors(modified_principal, rotation_matrix);
    mStressVector = this->RotatePrincipalStresses(modified_principal, rotation_matrix);
}

Vector MohrCoulombWithTensionCutOff::ReturnStressAtAxialZone(const Vector& rPrincipalTrialStressVector,
                                                             double TensileStrength) const
{
    Vector result = rPrincipalTrialStressVector;
    result(0)     = TensileStrength;
    result(2)     = TensileStrength - rPrincipalTrialStressVector(0) + rPrincipalTrialStressVector(2);
    return result;
}

Vector MohrCoulombWithTensionCutOff::ReturnStressAtCornerReturnZone(const Vector& rPrincipalTrialStressVector,
                                                                    const Vector& rCornerPoint) const
{
    Vector result = rPrincipalTrialStressVector;
    result(0)     = rCornerPoint(0) + rCornerPoint(1);
    result(2)     = rCornerPoint(0) - rCornerPoint(1);
    return result;
}

Vector MohrCoulombWithTensionCutOff::ReturnStressAtRegularFailureZone(const Vector& rPrincipalTrialStressVector,
                                                                      const CoulombYieldSurface& rCoulombYieldSurface,
                                                                      double FrictionAngle,
                                                                      double Cohesion) const
{
    Vector result = rPrincipalTrialStressVector;
    const auto flowDerivative = rCoulombYieldSurface.DerivateOfFlowFunction(rPrincipalTrialStressVector);
    const auto cof1   = (1.0 + std::sin(FrictionAngle)) / (1.0 - std::sin(FrictionAngle));
    const auto cof2   = 2.0 * Cohesion * std::cos(FrictionAngle) / (1.0 - std::sin(FrictionAngle));
    const auto det    = cof1 * flowDerivative(0) - flowDerivative(2);
    const auto lambda = (rPrincipalTrialStressVector(2) + cof2 - rPrincipalTrialStressVector(0) * cof1) / det;
    result(0) = lambda * flowDerivative(0) + rPrincipalTrialStressVector(0);
    result(2) = result(0) * cof1 - cof2;
    return result;
}

double MohrCoulombWithTensionCutOff::CalculateApex(double FrictionAngle, double Cohesion) const
{
    return Cohesion / std::tan(FrictionAngle);
}

Vector MohrCoulombWithTensionCutOff::CalculateCornerPoint(double FrictionAngle, double Cohesion, double TensileStrength) const
{
    Vector result = ZeroVector(2);
    if (TensileStrength < result(0)) {
        result(0) = this->CalculateApex(FrictionAngle, Cohesion);
        return result;
    }
    result(0) = (TensileStrength - Cohesion * std::cos(FrictionAngle)) / (1.0 - std::sin(FrictionAngle));
    result(1) = (Cohesion * std::cos(FrictionAngle) - TensileStrength * std::sin(FrictionAngle)) /
                (1.0 - std::sin(FrictionAngle));
    return result;
}

Vector MohrCoulombWithTensionCutOff::RotatePrincipalStresses(const Vector& rPrincipalStressVector,
                                                             const Matrix& rRotationMatrix) const
{
    Matrix principal_stress_matrix = GeoMechanicsMathUtilities::VectorToDiagonalMatrix(rPrincipalStressVector);
    Vector result = StressStrainUtilities::RotateStressMatrix(principal_stress_matrix, rRotationMatrix, mpConstitutiveDimension->GetStrainSize());
    return result;
}

void MohrCoulombWithTensionCutOff::CalculateTrialStressVector(const Vector& rStrainVector,
                                                              Vector&       rStressVector,
                                                              ConstitutiveLaw::Parameters& rValues) const
{
    const auto delta_strain_vector = rStrainVector - mStrainVectorFinalized;
    Matrix elastic_matrix    = this->CalculateElasticMatrix(rValues);
    rStressVector            = mStressVectorFinalized + prod(elastic_matrix, delta_strain_vector);
}

Matrix MohrCoulombWithTensionCutOff::CalculateElasticMatrix(ConstitutiveLaw::Parameters& rValues) const
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    return mpConstitutiveDimension->CalculateElasticMatrix(r_material_properties[YOUNG_MODULUS],
                                                           r_material_properties[POISSON_RATIO]);
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
}

void MohrCoulombWithTensionCutOff::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.load("StressVector", mStressVector);
    rSerializer.load("StressVectorFinalized", mStressVectorFinalized);
    rSerializer.load("StrainVectorFinalized", mStrainVectorFinalized);
}
} // Namespace Kratos