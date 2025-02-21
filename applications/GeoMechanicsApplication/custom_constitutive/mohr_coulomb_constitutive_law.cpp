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
#include "custom_constitutive/mohr_coulomb_constitutive_law.hpp"
#include "custom_constitutive/constitutive_law_dimension.h"
#include "custom_constitutive/coulomb_yield_function.hpp"
#include "custom_constitutive/tension_cutoff_function.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{

// MohrCoulombConstitutiveLaw::MohrCoulombConstitutiveLaw() = default;

MohrCoulombConstitutiveLaw::MohrCoulombConstitutiveLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mStressVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStressVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mDeltaStrainVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mStrainVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize()))
{
}

MohrCoulombConstitutiveLaw::~MohrCoulombConstitutiveLaw() = default;

ConstitutiveLaw::Pointer MohrCoulombConstitutiveLaw::Clone() const
{
    auto p_result = std::make_shared<MohrCoulombConstitutiveLaw>(mpConstitutiveDimension->Clone());
    p_result->mStressVector          = mStressVector;
    p_result->mStressVectorFinalized = mStressVectorFinalized;
    p_result->mDeltaStrainVector     = mDeltaStrainVector;
    p_result->mStrainVectorFinalized = mStrainVectorFinalized;
    return p_result;
}

int MohrCoulombConstitutiveLaw::Check(const Properties&   rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const ProcessInfo&  rCurrentProcessInfo) const
{
    // Verify Properties variables
    if (!rMaterialProperties.Has(GEO_COHESION) || rMaterialProperties[GEO_COHESION] <= 0.0)
        KRATOS_ERROR << "GEO_COHESION is not defined or has an invalid value for property: "
                     << rMaterialProperties.Id() << std::endl;

    if (!rMaterialProperties.Has(GEO_FRICTION_ANGLE) || rMaterialProperties[GEO_FRICTION_ANGLE] <= 0.0)
        KRATOS_ERROR << "GEO_FRICTION_ANGLE is not defined or has an invalid value for property: "
                     << rMaterialProperties.Id() << std::endl;

    if (!rMaterialProperties.Has(GEO_DILATION_ANGLE) || rMaterialProperties[GEO_DILATION_ANGLE] <= 0.0)
        KRATOS_ERROR << "GEO_DILATION_ANGLE is not defined or has an invalid value for property: "
                     << rMaterialProperties.Id() << std::endl;

    if (!rMaterialProperties.Has(GEO_TENSION_CUTOFF) || rMaterialProperties[GEO_TENSION_CUTOFF] <= 0.0)
        KRATOS_ERROR << "GEO_TENSION_CUTOFF is not defined or has an invalid value for property: "
                     << rMaterialProperties.Id() << std::endl;

    if (!rMaterialProperties.Has(YOUNG_MODULUS) || rMaterialProperties[YOUNG_MODULUS] <= 0.0)
        KRATOS_ERROR
            << "YOUNG_MODULUS has Key zero, is not defined or has an invalid value for property: "
            << rMaterialProperties.Id() << std::endl;

    if (!rMaterialProperties.Has(POISSON_RATIO) || rMaterialProperties[POISSON_RATIO] < 0.0)
        KRATOS_ERROR << "POISSON_RATIO is not defined or has an invalid value for property: "
                     << rMaterialProperties.Id() << std::endl;

    return 0;
}

double& MohrCoulombConstitutiveLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    if (rThisVariable == DAMAGE_VARIABLE || rThisVariable == STATE_VARIABLE)
        rValue = mStateVariable;
    return rValue;
}

void MohrCoulombConstitutiveLaw::SetValue(const Variable<double>& rThisVariable,
                                          const double&           rValue,
                                          const ProcessInfo&      rCurrentProcessInfo)
{
    if (rThisVariable == STATE_VARIABLE) mStateVariable = rValue;
}

void MohrCoulombConstitutiveLaw::CalculateMohrCoulomb(const Properties& rProp, Vector& rCauchyStressVector)
{
    const auto friction_angle = rProp[GEO_FRICTION_ANGLE] * Globals::Pi / 180.0;
    const auto cohesion       = rProp[GEO_COHESION];
    const auto dilation_angle = rProp[GEO_DILATION_ANGLE] * Globals::Pi / 180.0;
    const auto tension_cutoff = rProp[GEO_TENSION_CUTOFF];

    const auto coulombYieldFunction = CoulombYieldFunction(friction_angle, cohesion, dilation_angle);
    const auto tensionCutoffFunction = TensionCutoffFunction(tension_cutoff);

    ConstitutiveLaw::Parameters rValues ;
    rValues.SetMaterialProperties(rProp);
    Vector                       trailStressVector;
    Vector                       strainVector(mpConstitutiveDimension->GetStrainSize());
    rValues.SetStrainVector(strainVector);

    this->CalculatePK2Stress(strainVector, trailStressVector, rValues);
    Vector principalTrialStressVector = StressStrainUtilities::CalculatePrincipalStresses(trailStressVector);

    Matrix eigenVectorsMatrix = StressStrainUtilities::CalculatePrincipalEigenVectorsMatrix(trailStressVector);
    Matrix rotationMatrix = this->CalculateRotationMatrix(eigenVectorsMatrix);
    this->CheckRotationMatrix(rotationMatrix);

    double fme = coulombYieldFunction.CalculateYieldFunction(principalTrialStressVector);
    double fte = tensionCutoffFunction.CalculateYieldFunction(principalTrialStressVector);

    double trailStress = 0.5 * (principalTrialStressVector(0) + principalTrialStressVector(2));
    double trialShear  = 0.5 * (principalTrialStressVector(0) - principalTrialStressVector(2));

    // Elastic region
    if (fme <= 0.0 && fte<= 0.0) {
        mStressVector = this->ReturnStressAtElasticZone(trailStressVector);
        return;
    }

    // Zone of axial return
    double apex        = cohesion / std::tan(friction_angle);
    Vector cornerPoint = this->CalculateCornerPoint(friction_angle, cohesion, tension_cutoff, apex);

    if (tension_cutoff < apex && trialShear <= cornerPoint(1) && fte >= 0.0) {
        Vector modified_principal = this->ReturnStressAtAxialZone(principalTrialStressVector, tension_cutoff);
        mStressVector = this->RotatePrincipalStresses(modified_principal, rotationMatrix);
        return;
    }

    // Zone of tensile corner return
    double flow = (trailStress - cornerPoint(0)) - (trialShear - cornerPoint(1)) / std::sin(dilation_angle);
    if (flow <= 0.0 ) {    //&& trialShear > shearAtIntersection
        Vector modified_principal = this->ReturnStressAtCornerReturnZone(principalTrialStressVector, cornerPoint);
        mStressVector = this->RotatePrincipalStresses(modified_principal, rotationMatrix);
        return;
    }

    // Regular failure region
    if (fme > 0.0) { //&& flow > 0.0
        Vector modified_principal = this->ReturnStressAtRegularFailureZone(
            principalTrialStressVector, coulombYieldFunction, friction_angle, cohesion);
        mStressVector = this->RotatePrincipalStresses(modified_principal, rotationMatrix);
        return;
    }
}


// ================================================================================================
Vector MohrCoulombConstitutiveLaw::ReturnStressAtElasticZone(Vector& rTrailStressVector)
{
    return rTrailStressVector;
}

// ================================================================================================
Vector MohrCoulombConstitutiveLaw::ReturnStressAtAxialZone(Vector& rPrincipalTrialStressVector,
                                                           double  tension_cutoff)
{
    Vector result     = rPrincipalTrialStressVector;
    double trialShear = 0.5 * (rPrincipalTrialStressVector(0) - rPrincipalTrialStressVector(2));
    result(0) = tension_cutoff + 2.0 * trialShear;
    result(2) = tension_cutoff;
    return result;
}

// ================================================================================================
Vector MohrCoulombConstitutiveLaw::ReturnStressAtCornerReturnZone(Vector& rPrincipalTrialStressVector,
                                                                  Vector& rCornerPoint)
{
    Vector result = rPrincipalTrialStressVector;
    result(0) = rCornerPoint(0) + rCornerPoint(1);
    result(2) = rCornerPoint(0) - rCornerPoint(1);
    return result;
}

// ================================================================================================
Vector MohrCoulombConstitutiveLaw::ReturnStressAtRegularFailureZone(Vector& rPrincipalTrialStressVector,
                                                                    const CoulombYieldFunction& rCoulombYieldFunction,
                                                                    double  friction_angle,
                                                                    double  cohesion)
{
    Vector result = rPrincipalTrialStressVector;
    Vector flowDerivative = rCoulombYieldFunction.CalculateFlowFunctionDerivate(rPrincipalTrialStressVector);
    double cof1   = (1.0 + std::sin(friction_angle)) / (1.0 - std::sin(friction_angle));
    double cof2   = 2.0 * cohesion * std::cos(friction_angle) / (1.0 - std::sin(friction_angle));
    double det    = cof1 * flowDerivative(0) - flowDerivative(2);
    double lambda = rPrincipalTrialStressVector(2) + cof2 - rPrincipalTrialStressVector(0) * cof1;
    lambda /= det;
    result(0) = lambda * flowDerivative(0) + rPrincipalTrialStressVector(0);
    result(2) = result(0) * cof1 - cof2;
    return result;
}

// ================================================================================================
Vector MohrCoulombConstitutiveLaw::CalculateCornerPoint(double rFrictionAngle, double rCohesion, double rTensionCutoff, double rApex)
{
    Vector result = ZeroVector(2);
    result(1)     = 0.0;
    result(0)     = rApex;
    if (rTensionCutoff < rApex) {
        result(1) = (rCohesion * std::cos(rFrictionAngle) - rTensionCutoff * std::sin(rFrictionAngle)) /
                    (1.0 + std::sin(rFrictionAngle));
        result(0) = (rTensionCutoff + rCohesion * std::cos(rFrictionAngle)) / (1.0 + std::sin(rFrictionAngle));
    }
    return result;
}

// ================================================================================================
Vector MohrCoulombConstitutiveLaw::RotatePrincipalStresses(Vector& rPrincipalStressVector, Matrix& rRotationMatrix)
{
    Matrix principalStressMatrix = this->ConvertVectorToDiagonalMatrix(rPrincipalStressVector);
    Matrix temp         = prod(principalStressMatrix, trans(rRotationMatrix));
    Matrix stressMatrix = prod(rRotationMatrix, temp);
    return MathUtils<double>::StressTensorToVector(stressMatrix);
}

// ================================================================================================
void MohrCoulombConstitutiveLaw::CalculatePK2Stress(const Vector&                rStrainVector,
                                                    Vector&                      rStressVector,
                                                    ConstitutiveLaw::Parameters& rValues)
{
    mDeltaStrainVector = rValues.GetStrainVector() - mStrainVectorFinalized;

    Matrix C;
    this->CalculateElasticMatrix(C, rValues);

    // Incremental formulation
    auto trialStress = mStressVectorFinalized + prod(C, mDeltaStrainVector);

    rStressVector = trialStress;
}

// ================================================================================================
void MohrCoulombConstitutiveLaw::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const auto        E                     = r_material_properties[YOUNG_MODULUS];
    const auto        NU                    = r_material_properties[POISSON_RATIO];

    const double c0 = E / ((1.0 + NU) * (1.0 - 2.0 * NU));
    const double c1 = (1.0 - NU) * c0;
    const double c2 = c0 * NU;
    const double c3 = (0.5 - NU) * c0;

    C = mpConstitutiveDimension->FillConstitutiveMatrix(c1, c2, c3);
}

// ================================================================================================
void MohrCoulombConstitutiveLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mStrainVectorFinalized = rValues.GetStrainVector();
    mStressVectorFinalized = mStressVector;
}

// ================================================================================================
Vector MohrCoulombConstitutiveLaw::NormalizeVector(Vector& vector)
{ 
    double norm = MathUtils<>::Norm(vector);
    return vector / norm;
}

// ================================================================================================
Matrix MohrCoulombConstitutiveLaw::CalculateRotationMatrix(Matrix& eigenVectorsMatrix)
{
    Matrix result(eigenVectorsMatrix.size1(), eigenVectorsMatrix.size2());
    for (int i = 0; i < eigenVectorsMatrix.size1(); ++i) {
        Vector vec        = column(eigenVectorsMatrix, i);
        vec               = this->NormalizeVector(vec);
        column(result, i) = vec;
    }
    return result;
}

// ================================================================================================
void MohrCoulombConstitutiveLaw::CheckRotationMatrix(Matrix& rRotationMatrix)
{
    Matrix result = prod(rRotationMatrix, trans(rRotationMatrix));
    for (int i = 0; i < rRotationMatrix.size1(); ++i) {
        KRATOS_ERROR_IF_NOT(result(i, i) == 1.0) << "The rotation matrix is not orthogonal" << std::endl;
    }

    for (int i = 0; i < rRotationMatrix.size1(); ++i) {
        for (int j = 0; j < rRotationMatrix.size2(); ++j) {
            if (i != j)
                KRATOS_ERROR_IF_NOT(result(i, j) == 0.0)
                    << "The rotation matrix is not orthogonal" << std::endl;
        }
    }
}

// ================================================================================================
Matrix MohrCoulombConstitutiveLaw::ConvertVectorToDiagonalMatrix(Vector& rVector)
{
    Matrix result = ZeroMatrix(rVector.size(), rVector.size());
    for (int i = 0; i < rVector.size(); ++i) {
        result(i, i) = rVector(i);
    }
    return result;
}

} // Namespace Kratos