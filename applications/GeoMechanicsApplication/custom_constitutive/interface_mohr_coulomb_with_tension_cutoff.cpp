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
#include "custom_constitutive/interface_mohr_coulomb_with_tension_cutoff.h"
#include "custom_constitutive/constitutive_law_dimension.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/math_utilities.h"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

#include <cmath>

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

} // namespace

namespace Kratos
{
ConstitutiveLaw::Pointer InterfaceMohrCoulombWithTensionCutOff::Clone() const
{
    auto p_result             = std::make_shared<InterfaceMohrCoulombWithTensionCutOff>(*this);
    p_result->mTractionVector = mTractionVector;
    p_result->mTractionVectorFinalized             = mTractionVectorFinalized;
    p_result->mRelativeDisplacementVectorFinalized = mRelativeDisplacementVectorFinalized;
    p_result->mCoulombYieldSurface                 = mCoulombYieldSurface;
    p_result->mTensionCutOff                       = mTensionCutOff;
    return p_result;
}

Vector& InterfaceMohrCoulombWithTensionCutOff::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == CAUCHY_STRESS_VECTOR) {
        rValue = mTractionVector;
    } else {
        rValue = ConstitutiveLaw::GetValue(rThisVariable, rValue);
    }
    return rValue;
}

void InterfaceMohrCoulombWithTensionCutOff::SetValue(const Variable<Vector>& rVariable,
                                                     const Vector&           rValue,
                                                     const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == CAUCHY_STRESS_VECTOR) {
        mTractionVector = rValue;
    } else {
        KRATOS_ERROR << "Can't set value of " << rVariable.Name() << ": unsupported variable\n";
    }
}

SizeType InterfaceMohrCoulombWithTensionCutOff::WorkingSpaceDimension() { return 2; }

int InterfaceMohrCoulombWithTensionCutOff::Check(const Properties&   rMaterialProperties,
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
    CheckProperty(rMaterialProperties, INTERFACE_NORMAL_STIFFNESS);
    CheckProperty(rMaterialProperties, INTERFACE_SHEAR_STIFFNESS);
    return result;
}

ConstitutiveLaw::StressMeasure InterfaceMohrCoulombWithTensionCutOff::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

SizeType InterfaceMohrCoulombWithTensionCutOff::GetStrainSize() const
{
    return VOIGT_SIZE_2D_INTERFACE;
}

ConstitutiveLaw::StrainMeasure InterfaceMohrCoulombWithTensionCutOff::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

bool InterfaceMohrCoulombWithTensionCutOff::IsIncremental() { return true; }

bool InterfaceMohrCoulombWithTensionCutOff::RequiresInitializeMaterialResponse() { return true; }

void InterfaceMohrCoulombWithTensionCutOff::InitializeMaterial(const Properties& rMaterialProperties,
                                                               const Geometry<Node>&,
                                                               const Vector&)
{
    mCoulombWithTensionCutOffImpl = CoulombWithTensionCutOffImpl{
        MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE]), rMaterialProperties[GEO_COHESION],
        MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_DILATANCY_ANGLE]),
        rMaterialProperties[GEO_TENSILE_STRENGTH]};
    mCoulombYieldSurface =
        CoulombYieldSurface(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE]),
                            rMaterialProperties[GEO_COHESION],
                            MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_DILATANCY_ANGLE]));
    mTensionCutOff = TensionCutoff(rMaterialProperties[GEO_TENSILE_STRENGTH]);

    mRelativeDisplacementVectorFinalized =
        HasInitialState() ? GetInitialState().GetInitialStrainVector() : ZeroVector{GetStrainSize()};
    mTractionVectorFinalized =
        HasInitialState() ? GetInitialState().GetInitialStressVector() : ZeroVector{GetStrainSize()};
}

void InterfaceMohrCoulombWithTensionCutOff::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    if (!mIsModelInitialized) {
        mTractionVectorFinalized             = rValues.GetStressVector();
        mRelativeDisplacementVectorFinalized = rValues.GetStrainVector();
        mIsModelInitialized                  = true;
    }
}

void InterfaceMohrCoulombWithTensionCutOff::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rParameters)
{
    const auto& r_prop = rParameters.GetMaterialProperties();

    auto trial_sigma_tau = CalculateTrialTractionVector(
        rParameters.GetStrainVector(), r_prop[INTERFACE_NORMAL_STIFFNESS], r_prop[INTERFACE_SHEAR_STIFFNESS]);
    trial_sigma_tau[1] = std::abs(trial_sigma_tau[1]);

    if (!mCoulombWithTensionCutOffImpl.IsAdmissibleSigmaTau(trial_sigma_tau)) {
        trial_sigma_tau = ConstitutiveLawUtilities::ReturnMappingToCoulombWithTensionCutOff(
            r_prop, trial_sigma_tau, mCoulombYieldSurface, mTensionCutOff);
    }

    mTractionVector = trial_sigma_tau;

    rParameters.GetStressVector() = mTractionVector;
}

Vector InterfaceMohrCoulombWithTensionCutOff::CalculateTrialTractionVector(const Vector& rRelativeDisplacementVector,
                                                                           double NormalStiffness,
                                                                           double ShearStiffness) const
{
    return mTractionVectorFinalized + prod(MakeConstitutiveMatrix(NormalStiffness, ShearStiffness),
                                           rRelativeDisplacementVector - mRelativeDisplacementVectorFinalized);
}

void InterfaceMohrCoulombWithTensionCutOff::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mRelativeDisplacementVectorFinalized = rValues.GetStrainVector();
    mTractionVectorFinalized             = mTractionVector;
}

Matrix& InterfaceMohrCoulombWithTensionCutOff::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                              const Variable<Matrix>& rThisVariable,
                                                              Matrix&                 rValue)
{
    if (rThisVariable == CONSTITUTIVE_MATRIX) {
        const auto& r_properties = rParameterValues.GetMaterialProperties();
        rValue                   = MakeConstitutiveMatrix(r_properties[INTERFACE_NORMAL_STIFFNESS],
                                                          r_properties[INTERFACE_SHEAR_STIFFNESS]);
    } else {
        KRATOS_ERROR << "Can't calculate value of " << rThisVariable.Name() << ": unsupported variable\n";
    }

    return rValue;
}

Matrix InterfaceMohrCoulombWithTensionCutOff::MakeConstitutiveMatrix(double NormalStiffness, double ShearStiffness) const
{
    auto result  = Matrix{ZeroMatrix{GetStrainSize(), GetStrainSize()}};
    result(0, 0) = NormalStiffness;
    result(1, 1) = ShearStiffness;
    return result;
}

void InterfaceMohrCoulombWithTensionCutOff::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.save("TractionVector", mTractionVector);
    rSerializer.save("TractionVectorFinalized", mTractionVectorFinalized);
    rSerializer.save("RelativeDisplacementVectorFinalized", mRelativeDisplacementVectorFinalized);
    rSerializer.save("CoulombYieldSurface", mCoulombYieldSurface);
    rSerializer.save("TensionCutOff", mTensionCutOff);
    rSerializer.save("CoulombWithTensionCutOffImpl", mCoulombWithTensionCutOffImpl);
    rSerializer.save("IsModelInitialized", mIsModelInitialized);
}

void InterfaceMohrCoulombWithTensionCutOff::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.load("TractionVector", mTractionVector);
    rSerializer.load("TractionVectorFinalized", mTractionVectorFinalized);
    rSerializer.load("RelativeDisplacementVectorFinalized", mRelativeDisplacementVectorFinalized);
    rSerializer.load("CoulombYieldSurface", mCoulombYieldSurface);
    rSerializer.load("TensionCutOff", mTensionCutOff);
    rSerializer.load("CoulombWithTensionCutOffImpl", mCoulombWithTensionCutOffImpl);
    rSerializer.load("IsModelInitialized", mIsModelInitialized);
}

} // Namespace Kratos