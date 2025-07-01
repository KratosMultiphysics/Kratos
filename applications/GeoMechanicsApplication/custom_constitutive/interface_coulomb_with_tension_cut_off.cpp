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
#include "custom_constitutive/interface_coulomb_with_tension_cut_off.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_constants.h"
#include "geo_mechanics_application_variables.h"

#include <cmath>

namespace Kratos
{

ConstitutiveLaw::Pointer InterfaceCoulombWithTensionCutOff::Clone() const
{
    auto p_result                      = std::make_shared<InterfaceCoulombWithTensionCutOff>(*this);
    p_result->mTractionVector          = mTractionVector;
    p_result->mTractionVectorFinalized = mTractionVectorFinalized;
    p_result->mRelativeDisplacementVectorFinalized = mRelativeDisplacementVectorFinalized;
    p_result->mCoulombWithTensionCutOffImpl        = mCoulombWithTensionCutOffImpl;
    p_result->mIsModelInitialized                  = mIsModelInitialized;
    return p_result;
}

Vector& InterfaceCoulombWithTensionCutOff::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == CAUCHY_STRESS_VECTOR) {
        rValue = mTractionVector;
    } else {
        rValue = ConstitutiveLaw::GetValue(rVariable, rValue);
    }
    return rValue;
}

void InterfaceCoulombWithTensionCutOff::SetValue(const Variable<Vector>& rVariable,
                                                 const Vector&           rValue,
                                                 const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == CAUCHY_STRESS_VECTOR) {
        mTractionVector = rValue;
    } else {
        KRATOS_ERROR << "Can't set value of " << rVariable.Name() << ": unsupported variable\n";
    }
}

SizeType InterfaceCoulombWithTensionCutOff::WorkingSpaceDimension()
{
    // Note that this implementation assumes line interface elements. It needs to be modified when planar interface elements become available.
    return N_DIM_2D;
}

int InterfaceCoulombWithTensionCutOff::Check(const Properties&   rMaterialProperties,
                                             const GeometryType& rElementGeometry,
                                             const ProcessInfo&  rCurrentProcessInfo) const
{
    const auto result = ConstitutiveLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    ConstitutiveLawUtilities::CheckProperty(rMaterialProperties, GEO_COHESION);
    // Since `CheckProperty` doesn't accept a custom range yet, we have inlined a modified version
    // of `CheckProperty` to check the friction angle
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(GEO_FRICTION_ANGLE))
        << GEO_FRICTION_ANGLE.Name() << " is not defined for property " << rMaterialProperties.Id()
        << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[GEO_FRICTION_ANGLE] <= 0.0 || rMaterialProperties[GEO_FRICTION_ANGLE] >= 90.0)
        << "value of " << GEO_FRICTION_ANGLE.Name() << " for property " << rMaterialProperties.Id()
        << " is out of range: " << rMaterialProperties[GEO_FRICTION_ANGLE]
        << " is not in (0.0, 90.0)" << std::endl;
    ConstitutiveLawUtilities::CheckProperty(rMaterialProperties, GEO_DILATANCY_ANGLE,
                                            rMaterialProperties[GEO_FRICTION_ANGLE]);
    ConstitutiveLawUtilities::CheckProperty(
        rMaterialProperties, GEO_TENSILE_STRENGTH,
        rMaterialProperties[GEO_COHESION] /
            std::tan(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE])));
    ConstitutiveLawUtilities::CheckProperty(rMaterialProperties, INTERFACE_NORMAL_STIFFNESS);
    ConstitutiveLawUtilities::CheckProperty(rMaterialProperties, INTERFACE_SHEAR_STIFFNESS);
    return result;
}

ConstitutiveLaw::StressMeasure InterfaceCoulombWithTensionCutOff::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

SizeType InterfaceCoulombWithTensionCutOff::GetStrainSize() const
{
    // Note that this implementation assumes line interface elements. It needs to be modified when planar interface elements become available.
    return VOIGT_SIZE_2D_INTERFACE;
}

ConstitutiveLaw::StrainMeasure InterfaceCoulombWithTensionCutOff::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

bool InterfaceCoulombWithTensionCutOff::IsIncremental() { return true; }

bool InterfaceCoulombWithTensionCutOff::RequiresInitializeMaterialResponse() { return true; }

void InterfaceCoulombWithTensionCutOff::InitializeMaterial(const Properties& rMaterialProperties,
                                                           const Geometry<Node>&,
                                                           const Vector&)
{
    mCoulombWithTensionCutOffImpl = CoulombWithTensionCutOffImpl{
        MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE]), rMaterialProperties[GEO_COHESION],
        MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_DILATANCY_ANGLE]),
        rMaterialProperties[GEO_TENSILE_STRENGTH]};

    mRelativeDisplacementVectorFinalized =
        HasInitialState() ? GetInitialState().GetInitialStrainVector() : ZeroVector{GetStrainSize()};
    mTractionVectorFinalized =
        HasInitialState() ? GetInitialState().GetInitialStressVector() : ZeroVector{GetStrainSize()};
}

void InterfaceCoulombWithTensionCutOff::InitializeMaterialResponseCauchy(Parameters& rConstitutiveLawParameters)
{
    if (!mIsModelInitialized) {
        mTractionVectorFinalized             = rConstitutiveLawParameters.GetStressVector();
        mRelativeDisplacementVectorFinalized = rConstitutiveLawParameters.GetStrainVector();
        mIsModelInitialized                  = true;
    }
}

void InterfaceCoulombWithTensionCutOff::CalculateMaterialResponseCauchy(Parameters& rConstitutiveLawParameters)
{
    const auto& r_properties = rConstitutiveLawParameters.GetMaterialProperties();

    auto trial_sigma_tau = CalculateTrialTractionVector(rConstitutiveLawParameters.GetStrainVector(),
                                                        r_properties[INTERFACE_NORMAL_STIFFNESS],
                                                        r_properties[INTERFACE_SHEAR_STIFFNESS]);
    const auto negative = std::signbit(trial_sigma_tau[1]);
    trial_sigma_tau[1]  = std::abs(trial_sigma_tau[1]);

    if (!mCoulombWithTensionCutOffImpl.IsAdmissibleSigmaTau(trial_sigma_tau)) {
        trial_sigma_tau = mCoulombWithTensionCutOffImpl.DoReturnMapping(r_properties, trial_sigma_tau);
    }

    if (negative) trial_sigma_tau[1] *= -1.0;

    mTractionVector                              = trial_sigma_tau;
    rConstitutiveLawParameters.GetStressVector() = mTractionVector;
}

Vector InterfaceCoulombWithTensionCutOff::CalculateTrialTractionVector(const Vector& rRelativeDisplacementVector,
                                                                       double NormalStiffness,
                                                                       double ShearStiffness) const
{
    return mTractionVectorFinalized + prod(ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(
                                               NormalStiffness, ShearStiffness, GetStrainSize()),
                                           rRelativeDisplacementVector - mRelativeDisplacementVectorFinalized);
}

void InterfaceCoulombWithTensionCutOff::FinalizeMaterialResponseCauchy(Parameters& rConstitutiveLawParameters)
{
    mRelativeDisplacementVectorFinalized = rConstitutiveLawParameters.GetStrainVector();
    mTractionVectorFinalized             = mTractionVector;
}

Matrix& InterfaceCoulombWithTensionCutOff::CalculateValue(Parameters& rConstitutiveLawParameters,
                                                          const Variable<Matrix>& rVariable,
                                                          Matrix&                 rValue)
{
    if (rVariable == CONSTITUTIVE_MATRIX) {
        const auto& r_properties = rConstitutiveLawParameters.GetMaterialProperties();
        rValue                   = ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(
            r_properties[INTERFACE_NORMAL_STIFFNESS], r_properties[INTERFACE_SHEAR_STIFFNESS], GetStrainSize());
    } else {
        KRATOS_ERROR << "Can't calculate value of " << rVariable.Name() << ": unsupported variable\n";
    }

    return rValue;
}

void InterfaceCoulombWithTensionCutOff::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.save("TractionVector", mTractionVector);
    rSerializer.save("TractionVectorFinalized", mTractionVectorFinalized);
    rSerializer.save("RelativeDisplacementVectorFinalized", mRelativeDisplacementVectorFinalized);
    rSerializer.save("CoulombWithTensionCutOffImpl", mCoulombWithTensionCutOffImpl);
    rSerializer.save("IsModelInitialized", mIsModelInitialized);
}

void InterfaceCoulombWithTensionCutOff::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.load("TractionVector", mTractionVector);
    rSerializer.load("TractionVectorFinalized", mTractionVectorFinalized);
    rSerializer.load("RelativeDisplacementVectorFinalized", mRelativeDisplacementVectorFinalized);
    rSerializer.load("CoulombWithTensionCutOffImpl", mCoulombWithTensionCutOffImpl);
    rSerializer.load("IsModelInitialized", mIsModelInitialized);
}

} // Namespace Kratos