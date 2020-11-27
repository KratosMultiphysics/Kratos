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
#include "custom_utilities/fluid_element_utilities.h"
#include "fluid_dynamics_application_variables.h"

// Include base h
#include "qs_vms_residual_derivatives.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
int QSVMSResidualDerivatives<TDim, TNumNodes>::Check(
    const GeometryType& rGeometry,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    for (const auto& r_node : rGeometry) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RELAXED_ACCELERATION, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BODY_FORCE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADVPROJ, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DIVPROJ, r_node);
    }

    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod QSVMSResidualDerivatives<TDim, TNumNodes>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_2;
}

template <unsigned int TDim, unsigned int TNumNodes>
QSVMSResidualDerivatives<TDim, TNumNodes>::SecondDerivatives::SecondDerivatives(
    const Element& rElement,
    FluidConstitutiveLaw& rFluidConstitutiveLaw)
    : mrElement(rElement),
      mrFluidConstitutiveLaw(rFluidConstitutiveLaw)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void QSVMSResidualDerivatives<TDim, TNumNodes>::SecondDerivatives::Initialize(
    Matrix& rOutput,
    const ProcessInfo& rProcessInfo)
{
    const auto& r_geometry = mrElement.GetGeometry();

    const auto& properties = mrElement.GetProperties();
    mDensity = properties.GetValue(DENSITY);
    mDynamicViscosity = properties.GetValue(DYNAMIC_VISCOSITY);
    mElementSize = ElementSizeCalculator<TDim, TNumNodes>::MinimumElementSize(r_geometry);

    // setting up primal constitutive law
    InitializeConstitutiveLaw(mConstitutiveLawValues, mStrainRate, mShearStress, mC,
                              r_geometry, mrElement.GetProperties(), rProcessInfo);

    mDynamicTau = rProcessInfo[DYNAMIC_TAU];
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
        }
    }

    mBlockSize = rOutput.size2() / TNumNodes;
}

template <unsigned int TDim, unsigned int TNumNodes>
void QSVMSResidualDerivatives<TDim, TNumNodes>::SecondDerivatives::AddResidualDerivativeContributions(
    Matrix& rOutput,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX)
{
    using element_utilities = FluidElementUtilities<TNumNodes>;
    using derivative_utilities = QSVMSDerivativeUtilities<TDim>;

    const auto& r_geometry = mrElement.GetGeometry();

    array_1d<double, 3> velocity, mesh_velocity;

    // get gauss point evaluated values
    element_utilities::EvaluateInPoint(
        r_geometry, rN,
        std::tie(velocity, VELOCITY),
        std::tie(mesh_velocity, MESH_VELOCITY));

    const array_1d<double, 3> convective_velocity = velocity - mesh_velocity;
    const double convective_velocity_norm = norm_2(convective_velocity);

    BoundedVector<double, TNumNodes> convective_velocity_dot_dn_dx;
    element_utilities::Product(convective_velocity_dot_dn_dx, rdNdX, convective_velocity);

    double effective_viscosity;
    derivative_utilities::CalculateStrainRate(
        mStrainRate, mNodalVelocity, rdNdX);
    mConstitutiveLawValues.SetShapeFunctionsValues(rN);
    mrFluidConstitutiveLaw.CalculateMaterialResponseCauchy(mConstitutiveLawValues);
    mrFluidConstitutiveLaw.CalculateValue(
        mConstitutiveLawValues, EFFECTIVE_VISCOSITY, effective_viscosity);

    double tau_one, tau_two;
    CalculateTau(tau_one, tau_two, mElementSize, mDensity, effective_viscosity,
                 convective_velocity_norm, mDynamicTau, mDeltaTime);

    const double coeff_1 = W * mDensity;
    const double coeff_2 = coeff_1 * tau_one;

    // Note: Dof order is (u,v,[w,]p) for each node
    for (IndexType c = 0; c < TNumNodes; ++c) {
        const IndexType row = c * mBlockSize;
        for (IndexType a = 0; a < TNumNodes; ++a) {
            const IndexType col = a * mBlockSize;
            const double mass = coeff_1 * rN[a] * rN[c];
            for (IndexType d = 0; d < TDim; ++d) {
                double value = 0.0;

                value -= mass;
                value -= coeff_2 * mDensity * convective_velocity_dot_dn_dx[a] * rN[c];

                rOutput(row + d, col + d) += value;
                rOutput(row + d, col + TDim) -= coeff_2 * rdNdX(a, d) * rN[c];
            }
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void QSVMSResidualDerivatives<TDim, TNumNodes>::SecondDerivatives::Finalize(
    Matrix& rOutput,
    const ProcessInfo& rProcessInfo)
{
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
    InitializeConstitutiveLaw(mConstitutiveLawValues, mStrainRate, mShearStress, mC,
                              r_geometry, mrElement.GetProperties(), rProcessInfo);

    // setting up derivative constitutive law
    InitializeConstitutiveLaw(mConstitutiveLawValuesDerivative, mStrainRateDerivative,
                              mShearStressDerivative, mCDerivative, r_geometry,
                              mrElement.GetProperties(), rProcessInfo);

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
        r_geometry, rGaussPointShapeFunctions,
        std::tie(mPressure, PRESSURE),
        std::tie(mBodyForce, BODY_FORCE),
        std::tie(mVelocity, VELOCITY),
        std::tie(mMeshVelocity, MESH_VELOCITY),
        std::tie(mRelaxedAcceleration, RELAXED_ACCELERATION),
        std::tie(mMomentumProjection, ADVPROJ),
        std::tie(mMassProjection, DIVPROJ));

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

template <unsigned int TDim, unsigned int TNumNodes>
void QSVMSResidualDerivatives<TDim, TNumNodes>::InitializeConstitutiveLaw(
    ConstitutiveLaw::Parameters& rParameters,
    Vector& rStrainVector,
    Vector& rStressVector,
    Matrix& rConstitutiveMatrix,
    const GeometryType& rGeometry,
    const PropertiesType& rProperties,
    const ProcessInfo& rProcessInfo)
{
    rStrainVector.resize(TStrainSize);
    rStressVector.resize(TStrainSize);
    rConstitutiveMatrix.resize(TStrainSize, TStrainSize, false);

    rParameters = ConstitutiveLaw::Parameters(rGeometry, rProperties, rProcessInfo);

    auto& cl_options = rParameters.GetOptions();
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    rParameters.SetStrainVector(rStrainVector); // this is the input parameter
    rParameters.SetStressVector(rStressVector); // this is an ouput parameter
    rParameters.SetConstitutiveMatrix(rConstitutiveMatrix); // this is an ouput parameter
}

// template instantiations
template class QSVMSResidualDerivatives<2, 3>;
template class QSVMSResidualDerivatives<2, 4>;

template class QSVMSResidualDerivatives<3, 4>;
template class QSVMSResidualDerivatives<3, 8>;


} // namespace Kratos