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

void MohrCoulombConstitutiveLaw::CalculateMohrCoulomb(const Properties& rProp, Vector& rCautchyStressVector)
{
    Vector principalStress = StressStrainUtilities::CalculatePrincipalStresses(rCautchyStressVector);

    double friction_angle = rProp[GEO_FRICTION_ANGLE] * Globals::Pi / 180.0;
    double cohesion       = rProp[GEO_COHESION];
    double dilation_angle = rProp[GEO_DILATION_ANGLE] * Globals::Pi / 180.0;
    double tension_cutoff = rProp[GEO_TENSION_CUTOFF];

    auto const coulombYieldFunction = CoulombYieldFunction(friction_angle, cohesion, dilation_angle);
    auto const tensionCutoffFunction = TensionCutoffFunction(tension_cutoff);

    double fme = coulombYieldFunction.CalculateYieldFunction(principalStress);
    double fte = tensionCutoffFunction.CalculateYieldFunction(principalStress);

    int region_index = this->FindRegionIndex(fme, fte);
}

void MohrCoulombConstitutiveLaw::CalculatePK2Stress(const Vector&                rStrainVector,
                                                    Vector&                      rStressVector,
                                                    ConstitutiveLaw::Parameters& rValues)
{
    mDeltaStrainVector = rValues.GetStrainVector() - mStrainVectorFinalized;

    Matrix C;
    this->CalculateElasticMatrix(C, rValues);

    // Incremental formulation
    noalias(mStressVector) = mStressVectorFinalized + prod(C, mDeltaStrainVector);

    rStressVector = mStressVector;
}

void MohrCoulombConstitutiveLaw::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const auto        E                     = r_material_properties[YOUNG_MODULUS];
    const auto        NU                    = r_material_properties[POISSON_RATIO];

    const double c0 = E / ((1.0 + NU) * (1.0 - 2.0 * NU));
    const double c1 = (1.0 - NU) * c0;
    const double c2 = c0 * NU;
    const double c3 = (0.5 - NU) * c0;

    // C = mpConstitutiveDimension->FillConstitutiveMatrix(c1, c2, c3);
}

void MohrCoulombConstitutiveLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mStrainVectorFinalized = rValues.GetStrainVector();
    mStressVectorFinalized = mStressVector;
}

int MohrCoulombConstitutiveLaw::FindRegionIndex(double fme, double fte)
{
    if (fme <= 0.0 && fte <= 0.0) return 0;
    if (fme < 0.0 && fte <= 0.0) return 1;
    return 0;
}
} // Namespace Kratos