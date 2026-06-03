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
#include "custom_constitutive/interface_coulomb.h"
#include "custom_constitutive/constitutive_law_dimension.h"
#include "custom_constitutive/sigma_tau.hpp"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_constants.h"
#include "geo_mechanics_application_variables.h"

#include <cmath>

namespace Kratos
{

InterfaceCoulomb::InterfaceCoulomb(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mTractionVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mTractionVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mRelativeDisplacementVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize()))
{
}

ConstitutiveLaw::Pointer InterfaceCoulomb::Clone() const
{
    auto p_result = std::make_shared<InterfaceCoulomb>(mpConstitutiveDimension->Clone());
    p_result->mTractionVector                      = mTractionVector;
    p_result->mTractionVectorFinalized             = mTractionVectorFinalized;
    p_result->mRelativeDisplacementVectorFinalized = mRelativeDisplacementVectorFinalized;
    p_result->mCoulombImpl                         = mCoulombImpl;
    p_result->mIsModelInitialized                  = mIsModelInitialized;
    return p_result;
}

Vector& InterfaceCoulomb::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
        rValue = mTractionVector;
    } else {
        rValue = ConstitutiveLaw::GetValue(rVariable, rValue);
    }
    return rValue;
}

int& InterfaceCoulomb::GetValue(const Variable<int>& rVariable, int& rValue)
{
    if (rVariable == GEO_PLASTICITY_STATUS) {
        rValue = static_cast<int>(mCoulombImpl.GetPlasticityStatus());
    }
    return rValue;
}

void InterfaceCoulomb::SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
        mTractionVector = rValue;
    } else {
        KRATOS_ERROR << "Can't set value of " << rVariable.Name() << ": unsupported variable\n";
    }
}

SizeType InterfaceCoulomb::WorkingSpaceDimension()
{
    // Note that this implementation assumes line interface elements. It needs to be modified when planar interface elements become available.
    return N_DIM_2D;
}

int InterfaceCoulomb::Check(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const ProcessInfo&  rCurrentProcessInfo) const
{
    const auto result = ConstitutiveLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    const CheckProperties check_properties(rMaterialProperties, "property", CheckProperties::Bounds::AllInclusive);
    check_properties.Check(GEO_COHESION);
    constexpr auto max_value_angle = 90.0;
    check_properties.SingleUseBounds(CheckProperties::Bounds::AllExclusive).Check(GEO_FRICTION_ANGLE, 0.0, max_value_angle);
    check_properties.Check(GEO_DILATANCY_ANGLE, rMaterialProperties[GEO_FRICTION_ANGLE]);
    if (rMaterialProperties.Has(GEO_ENABLE_TENSION_CUT_OFF) && rMaterialProperties[GEO_ENABLE_TENSION_CUT_OFF]) {
        check_properties.Check(
            GEO_TENSILE_STRENGTH,
            rMaterialProperties[GEO_COHESION] /
                std::tan(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE])));
    }
    check_properties.Check(INTERFACE_NORMAL_STIFFNESS);
    check_properties.Check(INTERFACE_SHEAR_STIFFNESS);
    return result;
}

ConstitutiveLaw::StressMeasure InterfaceCoulomb::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

SizeType InterfaceCoulomb::GetStrainSize() const
{
    // Note that this implementation assumes line interface elements. It needs to be modified when planar interface elements become available.
    return VOIGT_SIZE_2D_INTERFACE;
}

ConstitutiveLaw::StrainMeasure InterfaceCoulomb::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

bool InterfaceCoulomb::IsIncremental() { return true; }

bool InterfaceCoulomb::RequiresInitializeMaterialResponse() { return true; }

void InterfaceCoulomb::InitializeMaterial(const Properties& rMaterialProperties, const Geometry<Node>&, const Vector&)
{
    mCoulombImpl = CoulombImpl{rMaterialProperties};

    mRelativeDisplacementVectorFinalized =
        HasInitialState() ? GetInitialState().GetInitialStrainVector() : ZeroVector{GetStrainSize()};
    mTractionVectorFinalized =
        HasInitialState() ? GetInitialState().GetInitialStressVector() : ZeroVector{GetStrainSize()};
}

void InterfaceCoulomb::InitializeMaterialResponseCauchy(Parameters& rConstitutiveLawParameters)
{
    if (!mIsModelInitialized) {
        mTractionVectorFinalized             = rConstitutiveLawParameters.GetStressVector();
        mRelativeDisplacementVectorFinalized = rConstitutiveLawParameters.GetStrainVector();
        mIsModelInitialized                  = true;
    }
}

void InterfaceCoulomb::CalculateMaterialResponseCauchy(Parameters& rConstitutiveLawParameters)
{
    if (!rConstitutiveLawParameters.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        return;
    }

    const auto& r_properties = rConstitutiveLawParameters.GetMaterialProperties();

    auto trial_sigma_tau = CalculateTrialTractionVector(rConstitutiveLawParameters.GetStrainVector(),
                                                        r_properties[INTERFACE_NORMAL_STIFFNESS],
                                                        r_properties[INTERFACE_SHEAR_STIFFNESS]);
    auto mapped_sigma_tau = trial_sigma_tau;

    const auto negative   = std::signbit(trial_sigma_tau.Tau());
    trial_sigma_tau.Tau() = std::abs(trial_sigma_tau.Tau());

    if (!mCoulombImpl.IsAdmissibleStressState(trial_sigma_tau)) {
        mapped_sigma_tau = mCoulombImpl.DoReturnMapping(
            trial_sigma_tau, mpConstitutiveDimension->CalculateElasticConstitutiveTensor(r_properties),
            Geo::PrincipalStresses::AveragingType::NO_AVERAGING);
        if (negative) mapped_sigma_tau.Tau() *= -1.0;
    }

    mTractionVector                              = mapped_sigma_tau.CopyTo<Vector>();
    rConstitutiveLawParameters.GetStressVector() = mTractionVector;
}

Geo::SigmaTau InterfaceCoulomb::CalculateTrialTractionVector(const Vector& rRelativeDisplacementVector,
                                                             double NormalStiffness,
                                                             double ShearStiffness) const
{
    constexpr auto number_of_normal_components = std::size_t{1};
    return Geo::SigmaTau{mTractionVectorFinalized +
                         prod(ConstitutiveLawUtilities::MakeInterfaceElasticConstitutiveTensor(
                                  NormalStiffness, ShearStiffness, GetStrainSize(), number_of_normal_components),
                              rRelativeDisplacementVector - mRelativeDisplacementVectorFinalized)};
}

void InterfaceCoulomb::FinalizeMaterialResponseCauchy(Parameters& rConstitutiveLawParameters)
{
    mRelativeDisplacementVectorFinalized = rConstitutiveLawParameters.GetStrainVector();
    mTractionVectorFinalized             = mTractionVector;
}

Matrix& InterfaceCoulomb::CalculateValue(Parameters&             rConstitutiveLawParameters,
                                         const Variable<Matrix>& rVariable,
                                         Matrix&                 rValue)
{
    if (rVariable == CONSTITUTIVE_MATRIX) {
        const auto&    r_properties = rConstitutiveLawParameters.GetMaterialProperties();
        constexpr auto number_of_normal_components = std::size_t{1};
        rValue = ConstitutiveLawUtilities::MakeInterfaceElasticConstitutiveTensor(
            r_properties[INTERFACE_NORMAL_STIFFNESS], r_properties[INTERFACE_SHEAR_STIFFNESS],
            GetStrainSize(), number_of_normal_components);
    } else {
        KRATOS_ERROR << "Can't calculate value of " << rVariable.Name() << ": unsupported variable\n";
    }

    return rValue;
}

void InterfaceCoulomb::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.save("ConstitutiveDimension", mpConstitutiveDimension);
    rSerializer.save("TractionVector", mTractionVector);
    rSerializer.save("TractionVectorFinalized", mTractionVectorFinalized);
    rSerializer.save("RelativeDisplacementVectorFinalized", mRelativeDisplacementVectorFinalized);
    rSerializer.save("CoulombImpl", mCoulombImpl);
    rSerializer.save("IsModelInitialized", mIsModelInitialized);
}

void InterfaceCoulomb::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.load("ConstitutiveDimension", mpConstitutiveDimension);
    rSerializer.load("TractionVector", mTractionVector);
    rSerializer.load("TractionVectorFinalized", mTractionVectorFinalized);
    rSerializer.load("RelativeDisplacementVectorFinalized", mRelativeDisplacementVectorFinalized);
    rSerializer.load("CoulombImpl", mCoulombImpl);
    rSerializer.load("IsModelInitialized", mIsModelInitialized);
}

} // Namespace Kratos