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
    if (rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
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
    if (rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
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

    const CheckProperties check_properties(rMaterialProperties, "property", CheckProperties::Bounds::AllInclusive);
    check_properties.Check(GEO_COHESION);
    constexpr auto max_value_angle = 90.0;
    check_properties.SingleUseBounds(CheckProperties::Bounds::AllExclusive).Check(GEO_FRICTION_ANGLE, 0.0, max_value_angle);
    check_properties.Check(GEO_DILATANCY_ANGLE, rMaterialProperties[GEO_FRICTION_ANGLE]);
    check_properties.Check(
        GEO_TENSILE_STRENGTH,
        rMaterialProperties[GEO_COHESION] /
            std::tan(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE])));
    check_properties.Check(INTERFACE_NORMAL_STIFFNESS);
    check_properties.Check(INTERFACE_SHEAR_STIFFNESS);
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
    mCoulombWithTensionCutOffImpl = CoulombWithTensionCutOffImpl{rMaterialProperties};

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
    auto mapped_sigma_tau = trial_sigma_tau;

    const auto negative   = std::signbit(trial_sigma_tau.Tau());
    trial_sigma_tau.Tau() = std::abs(trial_sigma_tau.Tau());

    if (!mCoulombWithTensionCutOffImpl.IsAdmissibleStressState(trial_sigma_tau)) {
        auto elastic_matrix = CalculateElasticContitativeMatrix(r_properties);
        mapped_sigma_tau    = mCoulombWithTensionCutOffImpl.DoReturnMapping(
            trial_sigma_tau, YieldSurface::YieldSurfaceAveragingType::NO_AVERAGING, elastic_matrix);
        if (negative) mapped_sigma_tau.Tau() *= -1.0;
    }

    mTractionVector                              = mapped_sigma_tau.CopyTo<Vector>();
    rConstitutiveLawParameters.GetStressVector() = mTractionVector;
}

Geo::SigmaTau InterfaceCoulombWithTensionCutOff::CalculateTrialTractionVector(const Vector& rRelativeDisplacementVector,
                                                                              double NormalStiffness,
                                                                              double ShearStiffness) const
{
    constexpr auto number_of_normal_components = std::size_t{1};
    return Geo::SigmaTau{mTractionVectorFinalized +
                         prod(ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(
                                  NormalStiffness, ShearStiffness, GetStrainSize(), number_of_normal_components),
                              rRelativeDisplacementVector - mRelativeDisplacementVectorFinalized)};
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
        const auto&    r_properties = rConstitutiveLawParameters.GetMaterialProperties();
        constexpr auto number_of_normal_components = std::size_t{1};
        rValue = ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(
            r_properties[INTERFACE_NORMAL_STIFFNESS], r_properties[INTERFACE_SHEAR_STIFFNESS],
            GetStrainSize(), number_of_normal_components);
    } else {
        KRATOS_ERROR << "Can't calculate value of " << rVariable.Name() << ": unsupported variable\n";
    }

    return rValue;
}

Matrix InterfaceCoulombWithTensionCutOff::CalculateElasticContitativeMatrix(const Properties& rMaterialProperties) const
{
    Matrix result = ZeroMatrix(2, 2);
    result(0, 0)  = rMaterialProperties[INTERFACE_NORMAL_STIFFNESS];
    result(1, 1)  = rMaterialProperties[INTERFACE_SHEAR_STIFFNESS];
    return result;
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