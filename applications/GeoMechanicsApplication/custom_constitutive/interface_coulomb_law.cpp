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
#include "custom_constitutive/constitutive_law_dimension.h"
#include "custom_constitutive/coulomb_impl.h"
#include "custom_constitutive/interface_coulomb_law.h"
#include "custom_constitutive/sigma_tau.hpp"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/math_utilities.hpp"
#include "geo_mechanics_application_constants.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{
InterfaceCoulombLaw::InterfaceCoulombLaw()                                          = default;
InterfaceCoulombLaw::~InterfaceCoulombLaw()                                         = default;
InterfaceCoulombLaw::InterfaceCoulombLaw(InterfaceCoulombLaw&&) noexcept            = default;
InterfaceCoulombLaw& InterfaceCoulombLaw::operator=(InterfaceCoulombLaw&&) noexcept = default;

InterfaceCoulombLaw::InterfaceCoulombLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mTractionVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mTractionVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mRelativeDisplacementVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mpCoulombImpl(std::make_unique<CoulombImpl>())
{
}

ConstitutiveLaw::Pointer InterfaceCoulombLaw::Clone() const
{
    auto p_result = std::make_shared<InterfaceCoulombLaw>(mpConstitutiveDimension->Clone());
    p_result->mTractionVector                      = mTractionVector;
    p_result->mTractionVectorFinalized             = mTractionVectorFinalized;
    p_result->mRelativeDisplacementVectorFinalized = mRelativeDisplacementVectorFinalized;
    p_result->mpCoulombImpl                        = mpCoulombImpl->Clone();
    p_result->mIsModelInitialized                  = mIsModelInitialized;
    p_result->mMaxRelativeOvershoot                = mMaxRelativeOvershoot;
    p_result->mMaxNumberOfSubSteps                 = mMaxNumberOfSubSteps;
    p_result->mCalculatedNumberOfSubSteps          = mCalculatedNumberOfSubSteps;
    return p_result;
}

Vector& InterfaceCoulombLaw::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
        rValue = mTractionVector;
    } else {
        rValue = ConstitutiveLaw::GetValue(rVariable, rValue);
    }
    return rValue;
}

int& InterfaceCoulombLaw::GetValue(const Variable<int>& rVariable, int& rValue)
{
    if (rVariable == GEO_PLASTICITY_STATUS) {
        rValue = static_cast<int>(mpCoulombImpl->GetPlasticityStatus());
    }
    return rValue;
}

void InterfaceCoulombLaw::SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
        mTractionVector = rValue;
    } else {
        KRATOS_ERROR << "Can't set value of " << rVariable.Name() << ": unsupported variable\n";
    }
}

SizeType InterfaceCoulombLaw::WorkingSpaceDimension()
{
    // Note that this implementation assumes line interface elements. It needs to be modified when planar interface elements become available.
    return N_DIM_2D;
}

int InterfaceCoulombLaw::Check(const Properties&   rMaterialProperties,
                               const GeometryType& rElementGeometry,
                               const ProcessInfo&  rCurrentProcessInfo) const
{
    const auto result = ConstitutiveLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    const CheckProperties check_properties(rMaterialProperties, "property", CheckProperties::Bounds::AllInclusive);
    check_properties.Check(GEO_COHESION);
    constexpr auto max_value_angle = 90.0;
    check_properties.SingleUseBounds(CheckProperties::Bounds::AllExclusive).Check(GEO_FRICTION_ANGLE, 0.0, max_value_angle);
    check_properties.Check(GEO_DILATANCY_ANGLE, rMaterialProperties[GEO_FRICTION_ANGLE]);
    if (ConstitutiveLawUtilities::WantTensionCutOff(rMaterialProperties)) {
        check_properties.Check(
            GEO_TENSILE_STRENGTH,
            rMaterialProperties[GEO_COHESION] /
                std::tan(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE])));
    }
    check_properties.Check(INTERFACE_NORMAL_STIFFNESS);
    check_properties.Check(INTERFACE_SHEAR_STIFFNESS);
    return result;
}

ConstitutiveLaw::StressMeasure InterfaceCoulombLaw::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

SizeType InterfaceCoulombLaw::GetStrainSize() const
{
    // Note that this implementation assumes line interface elements. It needs to be modified when planar interface elements become available.
    return VOIGT_SIZE_2D_INTERFACE;
}

ConstitutiveLaw::StrainMeasure InterfaceCoulombLaw::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

bool InterfaceCoulombLaw::IsIncremental() { return true; }

bool InterfaceCoulombLaw::RequiresInitializeMaterialResponse() { return true; }

void InterfaceCoulombLaw::InitializeMaterial(const Properties& rMaterialProperties, const Geometry<Node>&, const Vector&)
{
    mpCoulombImpl = std::make_unique<CoulombImpl>(rMaterialProperties);

    mRelativeDisplacementVectorFinalized =
        HasInitialState() ? GetInitialState().GetInitialStrainVector() : ZeroVector{GetStrainSize()};
    mTractionVectorFinalized =
        HasInitialState() ? GetInitialState().GetInitialStressVector() : ZeroVector{GetStrainSize()};

    if (rMaterialProperties.Has(GEO_MAX_RELATIVE_OVERSHOOT)) {
        mMaxRelativeOvershoot = rMaterialProperties[GEO_MAX_RELATIVE_OVERSHOOT];
    }
    if (rMaterialProperties.Has(GEO_MAX_NUMBER_OF_SUB_STEPS)) {
        mMaxNumberOfSubSteps = rMaterialProperties[GEO_MAX_NUMBER_OF_SUB_STEPS];
    }
}

void InterfaceCoulombLaw::InitializeMaterialResponseCauchy(Parameters& rConstitutiveLawParameters)
{
    if (!mIsModelInitialized) {
        mTractionVectorFinalized             = rConstitutiveLawParameters.GetStressVector();
        mRelativeDisplacementVectorFinalized = rConstitutiveLawParameters.GetStrainVector();
        mIsModelInitialized                  = true;
    }
}

void InterfaceCoulombLaw::CalculateMaterialResponseCauchy(Parameters& rConstitutiveLawParameters)
{
    if (!rConstitutiveLawParameters.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        return;
    }

    const auto& r_properties = rConstitutiveLawParameters.GetMaterialProperties();
    const auto elastic_matrix = mpConstitutiveDimension->CalculateElasticConstitutiveTensor(r_properties);
    const auto normal_stiffness = r_properties[INTERFACE_NORMAL_STIFFNESS];
    const auto shear_stiffness  = r_properties[INTERFACE_SHEAR_STIFFNESS];

    // Elastic predictor over the full increment starting from the last finalized state.
    const auto full_trial_sigma_tau = CalculateTrialTractionVector(
        rConstitutiveLawParameters.GetStrainVector(), normal_stiffness, shear_stiffness);

    // If the whole step stays elastic, there is no need to sub-step.
    auto aux_sigma_tau = full_trial_sigma_tau;
    aux_sigma_tau.Tau() = std::abs(aux_sigma_tau.Tau());
    if (mpCoulombImpl->IsAdmissibleStressState(aux_sigma_tau)) {
        mTractionVector                              = full_trial_sigma_tau.CopyTo<Vector>();
        rConstitutiveLawParameters.GetStressVector() = mTractionVector;
        return;
    }

    std::size_t number_of_sub_steps = 1;
    if (mMaxNumberOfSubSteps > 1) {
        // only calculate the number of sub-steps once per solution step for stability
        if (mCalculatedNumberOfSubSteps == 0) {
            mCalculatedNumberOfSubSteps = CalculateAdaptiveNumberOfSubSteps(aux_sigma_tau, elastic_matrix);
        }
        number_of_sub_steps = mCalculatedNumberOfSubSteps;
    }

    // Running committed state for the sub-stepping (start from last finalized state).
    Vector committed_traction     = mTractionVectorFinalized;
    Vector committed_displacement = mRelativeDisplacementVectorFinalized;

    const Vector total_displacement_increment =
        rConstitutiveLawParameters.GetStrainVector() - mRelativeDisplacementVectorFinalized;

    constexpr auto number_of_normal_components = std::size_t{1};
    const auto elastic_tensor = ConstitutiveLawUtilities::MakeInterfaceElasticConstitutiveTensor(
        normal_stiffness, shear_stiffness, GetStrainSize(), number_of_normal_components);

    for (std::size_t sub = 1; sub <= number_of_sub_steps; ++sub) {
        const Vector sub_displacement =
            mRelativeDisplacementVectorFinalized +
            (static_cast<double>(sub) / static_cast<double>(number_of_sub_steps)) * total_displacement_increment;

        auto trial_sigma_tau = Geo::SigmaTau{
            committed_traction + prod(elastic_tensor, sub_displacement - committed_displacement)};
        auto mapped_sigma_tau = trial_sigma_tau;

        const auto negative   = std::signbit(trial_sigma_tau.Tau());
        trial_sigma_tau.Tau() = std::abs(trial_sigma_tau.Tau());

        if (!mpCoulombImpl->IsAdmissibleStressState(trial_sigma_tau)) {
            mapped_sigma_tau = mpCoulombImpl->DoReturnMapping(
                trial_sigma_tau, elastic_matrix, Geo::PrincipalStresses::AveragingType::NO_AVERAGING);
            if (negative) mapped_sigma_tau.Tau() *= -1.0;
        }

        mTractionVector        = mapped_sigma_tau.CopyTo<Vector>();
        committed_traction     = mTractionVector;
        committed_displacement = sub_displacement;
    }

    rConstitutiveLawParameters.GetStressVector() = mTractionVector;
}

std::size_t InterfaceCoulombLaw::CalculateAdaptiveNumberOfSubSteps(const Geo::SigmaTau& rTrialTraction,
                                                                   const Matrix& rElasticMatrix)
{
    // make sure that kappa is not updated while calculating the number of required sub steps
    mpCoulombImpl->SaveKappaOfCoulombYieldSurface();
    const auto mapped_sigma_tau = mpCoulombImpl->DoReturnMapping(
        rTrialTraction, rElasticMatrix, Geo::PrincipalStresses::AveragingType::NO_AVERAGING);
    mpCoulombImpl->RestoreKappaOfCoulombYieldSurface();

    const Vector trial_values  = rTrialTraction.CopyTo<Vector>();
    const Vector mapped_values = mapped_sigma_tau.CopyTo<Vector>();

    const auto overshoot          = norm_2(trial_values - mapped_values);
    const auto stress_scale       = std::max(norm_2(trial_values), 1.0e-12);
    const auto relative_overshoot = overshoot / stress_scale;

    const auto number_of_sub_steps = static_cast<int>(std::ceil(relative_overshoot / mMaxRelativeOvershoot));

    return std::clamp(number_of_sub_steps, 1, mMaxNumberOfSubSteps);
}

Geo::SigmaTau InterfaceCoulombLaw::CalculateTrialTractionVector(const Vector& rRelativeDisplacementVector,
                                                                double NormalStiffness,
                                                                double ShearStiffness) const
{
    constexpr auto number_of_normal_components = std::size_t{1};
    return Geo::SigmaTau{mTractionVectorFinalized +
                         prod(ConstitutiveLawUtilities::MakeInterfaceElasticConstitutiveTensor(
                                  NormalStiffness, ShearStiffness, GetStrainSize(), number_of_normal_components),
                              rRelativeDisplacementVector - mRelativeDisplacementVectorFinalized)};
}

void InterfaceCoulombLaw::FinalizeMaterialResponseCauchy(Parameters& rConstitutiveLawParameters)
{
    mRelativeDisplacementVectorFinalized = rConstitutiveLawParameters.GetStrainVector();
    mTractionVectorFinalized             = mTractionVector;
    mCalculatedNumberOfSubSteps          = 0;
}

Matrix& InterfaceCoulombLaw::CalculateValue(Parameters&             rConstitutiveLawParameters,
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

void InterfaceCoulombLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.save("ConstitutiveDimension", mpConstitutiveDimension);
    rSerializer.save("TractionVector", mTractionVector);
    rSerializer.save("TractionVectorFinalized", mTractionVectorFinalized);
    rSerializer.save("RelativeDisplacementVectorFinalized", mRelativeDisplacementVectorFinalized);
    rSerializer.save("CoulombImpl", mpCoulombImpl);
    rSerializer.save("IsModelInitialized", mIsModelInitialized);
    rSerializer.save("MaxRelativeOvershoot", mMaxRelativeOvershoot);
    rSerializer.save("MaxNumberOfSubSteps", mMaxNumberOfSubSteps);
}

void InterfaceCoulombLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.load("ConstitutiveDimension", mpConstitutiveDimension);
    rSerializer.load("TractionVector", mTractionVector);
    rSerializer.load("TractionVectorFinalized", mTractionVectorFinalized);
    rSerializer.load("RelativeDisplacementVectorFinalized", mRelativeDisplacementVectorFinalized);
    rSerializer.load("CoulombImpl", mpCoulombImpl);
    rSerializer.load("IsModelInitialized", mIsModelInitialized);
    rSerializer.load("MaxRelativeOvershoot", mMaxRelativeOvershoot);
    rSerializer.load("MaxNumberOfSubSteps", mMaxNumberOfSubSteps);
}

} // Namespace Kratos
