//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"

// Project includes
#include "utilities/element_size_calculator.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

// Include base h
#include "qs_vms_residual_derivatives.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod QSVMSResidualDerivatives<TDim, TNumNodes>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_2;
}

template <unsigned int TDim, unsigned int TNumNodes>
QSVMSResidualDerivatives<TDim, TNumNodes>::Data::Data(
    const Element& rElement,
    FluidConstitutiveLaw& rFluidConstitutiveLaw)
    : mrElement(rElement),
      mrFluidConstitutiveLaw(rFluidConstitutiveLaw)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void QSVMSResidualDerivatives<TDim, TNumNodes>::Data::Initialize(const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    // this method gathers values which are constatnts for all gauss points.

    const auto& r_geometry = mrElement.GetGeometry();

    // get values from properties
    const auto& properties = mrElement.GetProperties();
    mDensity = properties.GetValue(DENSITY);
    mDynamicViscosity = properties.GetValue(DYNAMIC_VISCOSITY);

    // get values from process info
    mDynamicTau = rProcessInfo[DYNAMIC_TAU];
    mOSS_SWITCH = rProcessInfo[OSS_SWITCH];
    KRATOS_ERROR_IF(mOSS_SWITCH == 1)
        << "OSS Projection adjoints are not yet supported.\n";

    mDeltaTime = rProcessInfo[DELTA_TIME];
    KRATOS_ERROR_IF(mDeltaTime > 0.0)
        << "Adjoint is calculated in reverse time, therefore "
           "DELTA_TIME should be negative. [ DELTA_TIME = "
        << mDeltaTime << " ].\n";
    mDeltaTime *= -1.0;

    // filling nodal values
    for (IndexType a = 0; a < TNumNodes; ++a) {
        const auto& r_node = r_geometry[a];
        for (IndexType i = 0; i < TDim; ++i) {
            mNodalVelocity(a, i) = r_node.FastGetSolutionStepValue(VELOCITY)[i];
            mNodalMeshVelocity(a, i) =
                r_node.FastGetSolutionStepValue(MESH_VELOCITY)[i];
            mNodalEffectiveVelocity(a, i) =
                mNodalVelocity(a, i) - mNodalMeshVelocity(a, i);
        }

        mNodalPressure[a] = r_node.FastGetSolutionStepValue(PRESSURE);
    }

    // get other values
    mElementSize = ElementSizeCalculator<TDim, TNumNodes>::MinimumElementSize(r_geometry);

    // setting up primal constitutive law
    mStrainRate.resize(TStrainSize);
    mShearStress.resize(TStrainSize);
    mC.resize(TStrainSize, TStrainSize, false);

    mConstitutiveLawValues = ConstitutiveLaw::Parameters(
        r_geometry, mrElement.GetProperties(), rProcessInfo);

    auto& cl_options = mConstitutiveLawValues.GetOptions();
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    mConstitutiveLawValues.SetStrainVector(mStrainRate); // this is the input parameter
    mConstitutiveLawValues.SetStressVector(mShearStress); // this is an ouput parameter
    mConstitutiveLawValues.SetConstitutiveMatrix(mC); // this is an ouput parameter

    // setting up derivative constitutive law
    mStrainRateDerivative.resize(TStrainSize);
    mShearStressDerivative.resize(TStrainSize);
    mCDerivative.resize(TStrainSize, TStrainSize, false);

    mConstitutiveLawValuesDerivative = ConstitutiveLaw::Parameters(
        r_geometry, mrElement.GetProperties(), rProcessInfo);

    auto& cl_options_derivative = mConstitutiveLawValuesDerivative.GetOptions();
    cl_options_derivative.Set(ConstitutiveLaw::COMPUTE_STRESS);
    cl_options_derivative.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    mConstitutiveLawValuesDerivative.SetStrainVector(mStrainRateDerivative); // this is the input parameter
    mConstitutiveLawValuesDerivative.SetStressVector(mShearStressDerivative); // this is an ouput parameter
    mConstitutiveLawValuesDerivative.SetConstitutiveMatrix(mCDerivative); // this is an ouput parameter

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void QSVMSResidualDerivatives<TDim, TNumNodes>::Data::CalculateGaussPointData(
    const double GaussPointWeight,
    const Vector& rGaussPointShapeFunctions,
    const Matrix& rGaussPointShapeFunctionDerivatives)
{
    using element_utilities = FluidElementUtilities<TNumNodes>;
    using derivative_utilities = QSVMSDerivativeUtilities<TDim>;

    const auto& r_geometry = mrElement.GetGeometry();

    // get gauss point evaluated values
    element_utilities::EvaluateInPoint(
        r_geometry, rGaussPointShapeFunctions, std::tie(mPressure, PRESSURE),
        std::tie(mBodyForce, BODY_FORCE), std::tie(mVelocity, VELOCITY),
        std::tie(mMeshVelocity, MESH_VELOCITY),
        std::tie(mRelaxedAcceleration, RELAXED_ACCELERATION),
        std::tie(mMomentumProjection, ADVPROJ), std::tie(mMassProjection, DIVPROJ));

    mBodyForce *= mDensity;

    noalias(mConvectiveVelocity) = mVelocity - mMeshVelocity;
    mConvectiveVelocityNorm = norm_2(mConvectiveVelocity);
    element_utilities::Product(mConvectiveVelocityDotDnDx,
                               rGaussPointShapeFunctionDerivatives, mConvectiveVelocity);
    element_utilities::Product(mRelaxedAccelerationDotDnDx,
                               rGaussPointShapeFunctionDerivatives, mRelaxedAcceleration);
    element_utilities::Product(mBodyForceDotDnDx,
                               rGaussPointShapeFunctionDerivatives, mBodyForce);
    element_utilities::Product(mMomentumProjectionDotDnDx,
                               rGaussPointShapeFunctionDerivatives, mMomentumProjection);

    // compute constitutive law values
    // Ask Ruben: Why not (Velocity - MeshVelocity) in here?
    derivative_utilities::CalculateStrainRate(
        mStrainRate, mNodalVelocity, rGaussPointShapeFunctionDerivatives);
    mConstitutiveLawValues.SetShapeFunctionsValues(rGaussPointShapeFunctions);
    mConstitutiveLawValuesDerivative.SetShapeFunctionsValues(rGaussPointShapeFunctions);

    // ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
    // this is ok under the hypothesis that no history dependent behavior is employed
    mrFluidConstitutiveLaw.CalculateMaterialResponseCauchy(mConstitutiveLawValues);
    mrFluidConstitutiveLaw.CalculateValue(
        mConstitutiveLawValues, EFFECTIVE_VISCOSITY, mEffectiveViscosity);

    element_utilities::GetStrainMatrix(rGaussPointShapeFunctionDerivatives, mStrainMatrix);
    noalias(mViscousTermRHSContribution) = prod(trans(mStrainMatrix), mShearStress);

    CalculateTau(mTauOne, mTauTwo, mElementSize, mDensity, mEffectiveViscosity,
                 mConvectiveVelocityNorm, mDynamicTau, mDeltaTime);

    derivative_utilities::CalculateGradient(mPressureGradient, PRESSURE, r_geometry,
                                            rGaussPointShapeFunctionDerivatives);
    derivative_utilities::CalculateGradient(mVelocityGradient, VELOCITY, r_geometry,
                                            rGaussPointShapeFunctionDerivatives);
    derivative_utilities::CalculateGradient(mMeshVelocityGradient, MESH_VELOCITY, r_geometry,
                                            rGaussPointShapeFunctionDerivatives);
    noalias(mEffectiveVelocityGradient) = mVelocityGradient - mMeshVelocityGradient;

    element_utilities::Product(mPressureGradientDotDnDx,
                               rGaussPointShapeFunctionDerivatives, mPressureGradient);
    element_utilities::Product(mEffectiveVelocityDotVelocityGradient,
                               mVelocityGradient, mConvectiveVelocity);

    mVelocityDotNabla = 0.0;
    for (IndexType a = 0; a < TNumNodes; ++a) {
        for (IndexType i = 0; i < TDim; ++i) {
            mVelocityDotNabla +=
                rGaussPointShapeFunctionDerivatives(a, i) * mNodalVelocity(a, i);
        }
    }

    element_utilities::Product(mEffectiveVelocityDotVelocityGradientDotShapeGradient,
                               rGaussPointShapeFunctionDerivatives,
                               mEffectiveVelocityDotVelocityGradient);
}

template <unsigned int TDim, unsigned int TNumNodes>
void QSVMSResidualDerivatives<TDim, TNumNodes>::Data::Finalize(const ProcessInfo& rProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
double QSVMSResidualDerivatives<TDim, TNumNodes>::CalculateNormDerivative(
    const double ValueNorm,
    const Vector& Value,
    const Vector& ValueDerivative)
{
    if (ValueNorm > 0.0) {
        return inner_prod(Value, ValueDerivative) / ValueNorm;
    } else {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void QSVMSResidualDerivatives<TDim, TNumNodes>::CalculateTau(
    double& TauOne,
    double& TauTwo,
    const double ElementSize,
    const double Density,
    const double Viscosity,
    const double VelocityNorm,
    const double DynamicTau,
    const double DeltaTime)
{
    constexpr double c1 = 8.0;
    constexpr double c2 = 2.0;

    const double inv_tau =
        c1 * Viscosity / (ElementSize * ElementSize) +
        Density * (DynamicTau / DeltaTime + c2 * VelocityNorm / ElementSize);
    TauOne = 1.0 / inv_tau;
    TauTwo = Viscosity + c2 * Density * VelocityNorm * ElementSize / c1;
}

template <unsigned int TDim, unsigned int TNumNodes>
void QSVMSResidualDerivatives<TDim, TNumNodes>::CalculateTauDerivative(
    double& TauOneDerivative,
    double& TauTwoDerivative,
    const double TauOne,
    const double Density,
    const double DynamicTau,
    const double DeltaTime,
    const double ElementSize,
    const double ElementSizeDerivative,
    const double Viscosity,
    const double ViscosityDerivative,
    const double VelocityNorm,
    const double VelocityNormDerivative)
{
    constexpr double c1 = 8.0;
    constexpr double c2 = 2.0;

    const double h2 = std::pow(ElementSize, 2);
    const double h3 = std::pow(ElementSize, 3);

    double inv_tau_derivative = 0.0;
    inv_tau_derivative += c1 * ViscosityDerivative / h2;
    inv_tau_derivative += -2.0 * c1 * Viscosity * ElementSizeDerivative / h3;
    inv_tau_derivative += Density * c2 * VelocityNormDerivative / ElementSize;
    inv_tau_derivative += -1.0 * Density * c2 * VelocityNorm * ElementSizeDerivative / h2;
    TauOneDerivative = -1.0 * std::pow(TauOne, 2) * inv_tau_derivative;

    TauTwoDerivative = ViscosityDerivative;
    TauTwoDerivative += c2 * Density * VelocityNormDerivative * ElementSize / c1;
    TauTwoDerivative += c2 * Density * VelocityNorm * ElementSizeDerivative / c1;
}

// template instantiations

template class QSVMSResidualDerivatives<2, 3>;
template class QSVMSResidualDerivatives<2, 4>;

template class QSVMSResidualDerivatives<3, 4>;
template class QSVMSResidualDerivatives<3, 8>;


} // namespace Kratos