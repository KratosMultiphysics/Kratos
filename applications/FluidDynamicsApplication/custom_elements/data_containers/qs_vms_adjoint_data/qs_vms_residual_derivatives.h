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

#if !defined(KRATOS_QS_VMS_RESIDUAL_FIRST_DERIVATIVES_H)
#define KRATOS_QS_VMS_RESIDUAL_FIRST_DERIVATIVES_H

// System includes

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "utilities/time_discretization.h"

#include "utilities/element_size_calculator.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "custom_elements/data_containers/qs_vms_adjoint_data/qs_vms_derivative_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
class QSVMSResidualDerivatives
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    static constexpr IndexType TBlockSize = TDim + 1;

    /// Size of the strain and stress vectors (in Voigt notation) for the formulation
    constexpr static IndexType TStrainSize = (TDim - 1) * 3; // 3 in 2D, 6 in 3D

    constexpr static IndexType TElementLocalSize = TBlockSize * TNumNodes;

    ///@}
    ///@name Static Operations
    ///@{

    static GeometryData::IntegrationMethod GetIntegrationMethod()
    {
        return GeometryData::GI_GAUSS_2;
    }

    ///@}
    ///@name Class Declarations
    ///@{

    class Data;

    ///@}
    ///@name Classes
    ///@{

    template<class TDerivativesType, unsigned int TEquationOffset = 0>
    class VariableDerivatives
    {
    public:
        ///name@ Life Cycle
        ///@{

        VariableDerivatives(Data& rData) : mrData(rData) {}

        ///@}
        ///@name Operations
        ///@{

        void Initialize(
            const Matrix& rOutput,
            const ProcessInfo& rProcessInfo)
        {
            mBlockSize = rOutput.size2() / TNumNodes;
        }

        void CalculateResidualDerivative(
            BoundedVector<double, TElementLocalSize>& rResidualDerivative,
            const int NodeIndex,
            const int DirectionIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative)
        {
            using element_utilities = FluidElementUtilities<TNumNodes>;
            using derivative_utilities = QSVMSDerivativeUtilities<TDim>;

            rResidualDerivative.clear();

            constexpr double u_factor = TDerivativesType::VelocityDerivativeFactor;
            constexpr double p_factor = TDerivativesType::PressureDerivativeFactor;

            const auto& r_geometry = mrData.mrElement.GetGeometry();

            const TDerivativesType derivatives_type(
                NodeIndex, DirectionIndex, r_geometry, W, rN, rdNdX,
                WDerivative, DetJDerivative, rdNdXDerivative);

            const auto& velocity_derivative = derivatives_type.CalculateEffectiveVelocityDerivative(mrData.mVelocity);
            const auto& element_length_derivative = derivatives_type.CalculateElementLengthDerivative(mrData.mElementSize);
            derivatives_type.CalculateStrainRateDerivative(mrData.mStrainRateDerivative, mrData.mNodalVelocity);

            // compute viscous term derivatives
            auto p_fluid_constitutive_adjoint_law =
                mrData.mrFluidConstitutiveLaw.GetAdjointConstitutiveLaw();
            const double effective_viscosity_derivative =
                p_fluid_constitutive_adjoint_law->CalculateEffectiveViscosityDerivative(
                    mrData.mConstitutiveLawValuesDerivative, mrData.mConstitutiveLawValues,
                    NodeIndex, derivatives_type.GetDerivativeVariable());
            p_fluid_constitutive_adjoint_law->CalculateMaterialResponseCauchyDerivative(
                mrData.mConstitutiveLawValuesDerivative, mrData.mConstitutiveLawValues,
                NodeIndex, derivatives_type.GetDerivativeVariable(),
                mrData.mEffectiveViscosity, effective_viscosity_derivative);

            const double velocity_norm_derivative = CalculateNormDerivative(
                mrData.mConvectiveVelocityNorm, mrData.mConvectiveVelocity, velocity_derivative);

            double tau_one_derivative, tau_two_derivative;
            CalculateTauDerivative(
                tau_one_derivative, tau_two_derivative, mrData.mTauOne,
                mrData.mDensity, mrData.mDynamicTau, mrData.mDeltaTime,
                mrData.mElementSize, element_length_derivative,
                mrData.mEffectiveViscosity, effective_viscosity_derivative,
                mrData.mConvectiveVelocityNorm, velocity_norm_derivative);

            BoundedVector<double, TNumNodes> convective_velocity_derivative_dot_dn_dx;
            element_utilities::Product(
                convective_velocity_derivative_dot_dn_dx, rdNdX, velocity_derivative);
            BoundedVector<double, TNumNodes> relaxed_acceleration_dot_dn_dx_derivative;
            element_utilities::Product(relaxed_acceleration_dot_dn_dx_derivative, rdNdXDerivative, mrData.mRelaxedAcceleration);

            BoundedVector<double, TNumNodes> convective_velocity_dot_dn_dx_derivative;
            element_utilities::Product(
                convective_velocity_dot_dn_dx_derivative, rdNdXDerivative,
                mrData.mConvectiveVelocity);

            array_1d<double, 3> pressure_gradient_derivative;
            derivative_utilities::CalculateGradient(pressure_gradient_derivative, PRESSURE, r_geometry, rdNdXDerivative);
            for (IndexType l = 0; l < TDim; ++l) {
                pressure_gradient_derivative[l] += rdNdX(NodeIndex, l) * p_factor;
            }

            BoundedMatrix<double, TDim, TDim> velocity_gradient_derivative;
            element_utilities::Product(velocity_gradient_derivative,
                                       trans(mrData.mNodalVelocity), rdNdXDerivative);
            row(velocity_gradient_derivative, DirectionIndex) += row(rdNdX, NodeIndex) * u_factor;

            array_1d<double, 3> effective_velocity_derivative_dot_velocity_gradient;
            element_utilities::Product(effective_velocity_derivative_dot_velocity_gradient,
                                       mrData.mVelocityGradient, velocity_derivative);

            array_1d<double, 3> effective_velocity_dot_velocity_gradient_derivative;
            element_utilities::Product(
                effective_velocity_dot_velocity_gradient_derivative,
                velocity_gradient_derivative, mrData.mConvectiveVelocity);

            double velocity_dot_nabla_derivative = 0.0;
            for (IndexType a = 0; a < TNumNodes; ++a) {
                for (IndexType i = 0; i < TDim; ++i) {
                    velocity_dot_nabla_derivative += mrData.mNodalVelocity(a, i) * rdNdXDerivative(a, i);
                }
            }

            BoundedVector<double, TNumNodes> effective_velocity_derivative_dot_velocity_gradient_dot_shape_gradient;
            element_utilities::Product(
                effective_velocity_derivative_dot_velocity_gradient_dot_shape_gradient,
                rdNdX, effective_velocity_derivative_dot_velocity_gradient);

            BoundedVector<double, TNumNodes> effective_velocity_dot_velocity_gradient_derivative_dot_shape_gradient;
            element_utilities::Product(
                effective_velocity_dot_velocity_gradient_derivative_dot_shape_gradient,
                rdNdX, effective_velocity_dot_velocity_gradient_derivative);

            BoundedVector<double, TNumNodes> effective_velocity_dot_velocity_gradient_dot_shape_gradient_derivative;
            element_utilities::Product(
                effective_velocity_dot_velocity_gradient_dot_shape_gradient_derivative,
                rdNdXDerivative, mrData.mEffectiveVelocityDotVelocityGradient);

            BoundedVector<double, TNumNodes> pressure_gradient_derivative_dot_shape_gradient;
            element_utilities::Product(pressure_gradient_derivative_dot_shape_gradient,
                                       rdNdX, pressure_gradient_derivative);

            BoundedVector<double, TNumNodes> pressure_gradient_dot_shape_gradient_derivative;
            element_utilities::Product(pressure_gradient_dot_shape_gradient_derivative,
                                       rdNdXDerivative, mrData.mPressureGradient);

            // TODO: Needs implementation for OSS projections
            const array_1d<double, 3> momentum_projection_derivative = ZeroVector(3);
            const double mass_projection_derivative = 0.0;

            for (IndexType a = 0; a < TNumNodes; ++a) {
                const IndexType row = a * mBlockSize + TEquationOffset;

                double forcing_derivative = 0.0;
                for (IndexType i = 0; i < TDim; ++i) {

                    double value = 0.0;

                    // Adding RHS derivative terms
                    value += WDerivative * rN[a] * mrData.mBodyForce[i]; // v*BodyForce

                    value += mrData.mDensity * WDerivative * mrData.mTauOne * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mBodyForce[i];
                    value += mrData.mDensity * W * tau_one_derivative * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mBodyForce[i];
                    value += mrData.mDensity * W * mrData.mTauOne * convective_velocity_derivative_dot_dn_dx[a] * mrData.mBodyForce[i];
                    value += mrData.mDensity * W * mrData.mTauOne * convective_velocity_dot_dn_dx_derivative[a] * mrData.mBodyForce[i];

                    value -= mrData.mDensity * WDerivative * mrData.mTauOne * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mMomentumProjection[i];
                    value -= mrData.mDensity * W * tau_one_derivative * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mMomentumProjection[i];
                    value -= mrData.mDensity * W * mrData.mTauOne * convective_velocity_derivative_dot_dn_dx[a] * mrData.mMomentumProjection[i];
                    value -= mrData.mDensity * W * mrData.mTauOne * convective_velocity_dot_dn_dx_derivative[a] * mrData.mMomentumProjection[i];
                    value -= mrData.mDensity * W * mrData.mTauOne * mrData.mConvectiveVelocityDotDnDx[a] * momentum_projection_derivative[i];

                    value -= WDerivative * mrData.mTauTwo * rdNdX(a, i) * mrData.mMassProjection;
                    value -= W * tau_two_derivative * rdNdX(a, i) * mrData.mMassProjection;
                    value -= W * mrData.mTauTwo * rdNdXDerivative(a, i) * mrData.mMassProjection;
                    value -= W * mrData.mTauTwo * rdNdX(a, i) * mass_projection_derivative;

                    forcing_derivative += rdNdXDerivative(a, i) * mrData.mBodyForce[i];
                    forcing_derivative -= rdNdXDerivative(a, i) * mrData.mMomentumProjection[i];
                    forcing_derivative -= rdNdX(a, i) * momentum_projection_derivative[i];

                    // Adding LHS derivative terms
                    value -= WDerivative * mrData.mDensity * rN[a] * mrData.mEffectiveVelocityDotVelocityGradient[i];
                    value -= W * mrData.mDensity * rN[a] * effective_velocity_derivative_dot_velocity_gradient[i];
                    value -= W * mrData.mDensity * rN[a] * effective_velocity_dot_velocity_gradient_derivative[i];

                    value -= WDerivative * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mTauOne * mrData.mDensity * mrData.mEffectiveVelocityDotVelocityGradient[i];
                    value -= W * mrData.mDensity * convective_velocity_derivative_dot_dn_dx[a] * mrData.mTauOne * mrData.mDensity * mrData.mEffectiveVelocityDotVelocityGradient[i];
                    value -= W * mrData.mDensity * convective_velocity_dot_dn_dx_derivative[a] * mrData.mTauOne * mrData.mDensity * mrData.mEffectiveVelocityDotVelocityGradient[i];
                    value -= W * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * tau_one_derivative * mrData.mDensity * mrData.mEffectiveVelocityDotVelocityGradient[i];
                    value -= W * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mTauOne * mrData.mDensity * effective_velocity_derivative_dot_velocity_gradient[i];
                    value -= W * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mTauOne * mrData.mDensity * effective_velocity_dot_velocity_gradient_derivative[i];

                    value -= WDerivative * mrData.mTauOne * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mPressureGradient[i];
                    value -= W * tau_one_derivative * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mPressureGradient[i];
                    value -= W * mrData.mTauOne * mrData.mDensity * convective_velocity_derivative_dot_dn_dx[a] * mrData.mPressureGradient[i];
                    value -= W * mrData.mTauOne * mrData.mDensity * convective_velocity_dot_dn_dx_derivative[a] * mrData.mPressureGradient[i];
                    value -= W * mrData.mTauOne * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * pressure_gradient_derivative[i];

                    value += WDerivative * rdNdX(a, i) * mrData.mPressure;
                    value += W * rdNdXDerivative(a, i) * mrData.mPressure;
                    value += W * rdNdX(a, i) * rN[NodeIndex] * p_factor;

                    value -= WDerivative * mrData.mTauTwo * rdNdX(a, i) * mrData.mVelocityDotNabla;
                    value -= W * tau_two_derivative * rdNdX(a, i) * mrData.mVelocityDotNabla;
                    value -= W * mrData.mTauTwo * rdNdXDerivative(a, i) * mrData.mVelocityDotNabla;
                    value -= W * mrData.mTauTwo * rdNdX(a, i) * velocity_dot_nabla_derivative;
                    value -= W * mrData.mTauTwo * rdNdX(a, i) * rdNdX(NodeIndex, DirectionIndex) * u_factor;

                    // Adding Mass term derivatives
                    value -= WDerivative * mrData.mDensity * rN[a] * mrData.mRelaxedAcceleration[i];

                    value -= WDerivative * mrData.mTauOne * mrData.mDensity * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mRelaxedAcceleration[i];
                    value -= W * tau_one_derivative * mrData.mDensity * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mRelaxedAcceleration[i];
                    value -= W * mrData.mTauOne * mrData.mDensity * mrData.mDensity * convective_velocity_derivative_dot_dn_dx[a] * mrData.mRelaxedAcceleration[i];
                    value -= W * mrData.mTauOne * mrData.mDensity * mrData.mDensity * convective_velocity_dot_dn_dx_derivative[a] * mrData.mRelaxedAcceleration[i];

                    rResidualDerivative[row + i] += value;
                }

                double value = 0.0;

                const double forcing = mrData.mBodyForceDotDnDx[a] - mrData.mMomentumProjectionDotDnDx[a];
                value += WDerivative * mrData.mTauOne * forcing;
                value += W * tau_one_derivative * forcing;
                value += W * mrData.mTauOne * forcing_derivative;

                value -= WDerivative * mrData.mTauOne * mrData.mDensity * mrData.mEffectiveVelocityDotVelocityGradientDotShapeGradient[a];
                value -= W * tau_one_derivative * mrData.mDensity * mrData.mEffectiveVelocityDotVelocityGradientDotShapeGradient[a];
                value -= W * mrData.mTauOne * mrData.mDensity * effective_velocity_derivative_dot_velocity_gradient_dot_shape_gradient[a];
                value -= W * mrData.mTauOne * mrData.mDensity * effective_velocity_dot_velocity_gradient_derivative_dot_shape_gradient[a];
                value -= W * mrData.mTauOne * mrData.mDensity * effective_velocity_dot_velocity_gradient_dot_shape_gradient_derivative[a];

                value -= WDerivative * rN[a] * mrData.mVelocityDotNabla;
                value -= W * rN[a] * velocity_dot_nabla_derivative;
                value -= W * rN[a] * rdNdX(NodeIndex, DirectionIndex) * u_factor;

                value -= WDerivative * mrData.mTauOne * mrData.mPressureGradientDotDnDx[a];
                value -= W * tau_one_derivative * mrData.mPressureGradientDotDnDx[a];
                value -= W * mrData.mTauOne * pressure_gradient_derivative_dot_shape_gradient[a];
                value -= W * mrData.mTauOne * pressure_gradient_dot_shape_gradient_derivative[a];

                // Adding mass term derivatives
                value -=  WDerivative * mrData.mTauOne * mrData.mDensity * mrData.mRelaxedAccelerationDotDnDx[a];
                value -=  W * tau_one_derivative * mrData.mDensity * mrData.mRelaxedAccelerationDotDnDx[a];
                value -=  W * mrData.mTauOne * mrData.mDensity * relaxed_acceleration_dot_dn_dx_derivative[a];

                rResidualDerivative[row + TDim] += value;
            }

            this->AddViscousDerivative(rResidualDerivative, NodeIndex,
                                       DirectionIndex, W, rN, rdNdX, WDerivative,
                                       DetJDerivative, rdNdXDerivative);
        }

        void Finalize(
            const Matrix& rOutput,
            const ProcessInfo& rProcessInfo)
        {
        }

        ///@}

    private:
        ///@name Private Members
        ///@{

        Data& mrData;
        IndexType mBlockSize;

        ///@}
        ///@name Private Operations
        ///@{

        void AddViscousDerivative(
            BoundedVector<double, TElementLocalSize>& rResidualDerivative,
            const int NodeIndex,
            const int DirectionIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative) const
        {
            BoundedMatrix<double, TStrainSize, TElementLocalSize> strain_matrix_derivative;
            FluidElementUtilities<TNumNodes>::GetStrainMatrix(rdNdXDerivative, strain_matrix_derivative);

            BoundedVector<double, TElementLocalSize> rhs_contribution_derivative =
                prod(trans(mrData.mStrainMatrix), mrData.mShearStressDerivative) +
                prod(trans(strain_matrix_derivative), mrData.mShearStress);

            for (IndexType a = 0; a < TNumNodes; ++a) {
                const IndexType row = a * mBlockSize + TEquationOffset;
                const IndexType local_row = a * TBlockSize;

                for (IndexType i = 0; i < TDim; ++i) {
                    double value = 0.0;

                    value -= WDerivative * mrData.mViscousTermRHSContribution[local_row + i];
                    value -= W * rhs_contribution_derivative[local_row + i];

                    rResidualDerivative[row + i] += value;
                }

                double value = 0.0;

                value -= WDerivative * mrData.mViscousTermRHSContribution[local_row + TDim];
                value -= W * rhs_contribution_derivative[local_row + TDim];

                rResidualDerivative[row + TDim] += value;
            }
        }

        ///@}
    };

    class Data
    {
    public:
        ///@name Life Cycle
        ///@{

        Data(
            const Element& rElement,
            FluidConstitutiveLaw& rFluidConstitutiveLaw)
            : mrElement(rElement),
              mrFluidConstitutiveLaw(rFluidConstitutiveLaw)
        {
        }

        ///@}
        ///@name Operations
        ///@{

        void Initialize(const ProcessInfo& rProcessInfo)
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
                    mNodalMeshVelocity(a, i) = r_node.FastGetSolutionStepValue(MESH_VELOCITY)[i];
                    mNodalEffectiveVelocity(a, i) = mNodalVelocity(a, i) - mNodalMeshVelocity(a, i);
                }

                mNodalPressure[a] = r_node.FastGetSolutionStepValue(PRESSURE);
            }

            // get other values
            mElementSize = ElementSizeCalculator<TDim,TNumNodes>::MinimumElementSize(r_geometry);

            // setting up primal constitutive law
            mStrainRate.resize(TStrainSize);
            mShearStress.resize(TStrainSize);
            mC.resize(TStrainSize, TStrainSize, false);

            mConstitutiveLawValues = ConstitutiveLaw::Parameters(r_geometry, mrElement.GetProperties(), rProcessInfo);

            auto& cl_options = mConstitutiveLawValues.GetOptions();
            cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
            cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

            mConstitutiveLawValues.SetStrainVector(mStrainRate);   //this is the input parameter
            mConstitutiveLawValues.SetStressVector(mShearStress);  //this is an ouput parameter
            mConstitutiveLawValues.SetConstitutiveMatrix(mC);      //this is an ouput parameter

            // setting up derivative constitutive law
            mStrainRateDerivative.resize(TStrainSize);
            mShearStressDerivative.resize(TStrainSize);
            mCDerivative.resize(TStrainSize, TStrainSize, false);

            mConstitutiveLawValuesDerivative = ConstitutiveLaw::Parameters(r_geometry, mrElement.GetProperties(), rProcessInfo);

            auto& cl_options_derivative = mConstitutiveLawValuesDerivative.GetOptions();
            cl_options_derivative.Set(ConstitutiveLaw::COMPUTE_STRESS);
            cl_options_derivative.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

            mConstitutiveLawValuesDerivative.SetStrainVector(mStrainRateDerivative);   //this is the input parameter
            mConstitutiveLawValuesDerivative.SetStressVector(mShearStressDerivative);  //this is an ouput parameter
            mConstitutiveLawValuesDerivative.SetConstitutiveMatrix(mCDerivative);      //this is an ouput parameter

            KRATOS_CATCH("");
        }

        void CalculateGaussPointData(
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
            element_utilities::Product(mConvectiveVelocityDotDnDx, rGaussPointShapeFunctionDerivatives, mConvectiveVelocity);
            element_utilities::Product(mRelaxedAccelerationDotDnDx, rGaussPointShapeFunctionDerivatives, mRelaxedAcceleration);
            element_utilities::Product(mBodyForceDotDnDx, rGaussPointShapeFunctionDerivatives, mBodyForce);
            element_utilities::Product(mMomentumProjectionDotDnDx, rGaussPointShapeFunctionDerivatives, mMomentumProjection);

            // compute constitutive law values
            // Ask Ruben: Why not (Velocity - MeshVelocity) in here?
            derivative_utilities::CalculateStrainRate(mStrainRate, mNodalVelocity, rGaussPointShapeFunctionDerivatives);
            mConstitutiveLawValues.SetShapeFunctionsValues(rGaussPointShapeFunctions);
            mConstitutiveLawValuesDerivative.SetShapeFunctionsValues(rGaussPointShapeFunctions);

            //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
            //this is ok under the hypothesis that no history dependent behavior is employed
            mrFluidConstitutiveLaw.CalculateMaterialResponseCauchy(mConstitutiveLawValues);
            mrFluidConstitutiveLaw.CalculateValue(mConstitutiveLawValues, EFFECTIVE_VISCOSITY, mEffectiveViscosity);

            element_utilities::GetStrainMatrix(rGaussPointShapeFunctionDerivatives, mStrainMatrix);
            noalias(mViscousTermRHSContribution) = prod(trans(mStrainMatrix), mShearStress) * (-1.0);

            CalculateTau(mTauOne, mTauTwo, mElementSize, mDensity, mEffectiveViscosity,
                         mConvectiveVelocityNorm, mDynamicTau, mDeltaTime);

            derivative_utilities::CalculateGradient(mPressureGradient, PRESSURE, r_geometry, rGaussPointShapeFunctionDerivatives);
            derivative_utilities::CalculateGradient(mVelocityGradient, VELOCITY, r_geometry, rGaussPointShapeFunctionDerivatives);
            derivative_utilities::CalculateGradient(mMeshVelocityGradient, MESH_VELOCITY, r_geometry, rGaussPointShapeFunctionDerivatives);
            noalias(mEffectiveVelocityGradient) = mVelocityGradient - mMeshVelocityGradient;

            element_utilities::Product(mPressureGradientDotDnDx, rGaussPointShapeFunctionDerivatives, mPressureGradient);
            element_utilities::Product(mEffectiveVelocityDotVelocityGradient,
                                       mVelocityGradient, mConvectiveVelocity);

            mVelocityDotNabla = 0.0;
            for (IndexType a = 0; a < TNumNodes; ++a) {
                for (IndexType i = 0; i < TDim; ++i) {
                    mVelocityDotNabla += rGaussPointShapeFunctionDerivatives(a, i) * mNodalVelocity(a, i);
                }
            }

            element_utilities::Product(mEffectiveVelocityDotVelocityGradientDotShapeGradient,
                                       rGaussPointShapeFunctionDerivatives,
                                       mEffectiveVelocityDotVelocityGradient);
        }

        void Finalize(const ProcessInfo& rProcessInfo)
        {
        }

        ///@}
    private:
        ///@name Private Members
        ///@{

        const Element& mrElement;
        FluidConstitutiveLaw& mrFluidConstitutiveLaw;

        // Primal data
        int mOSS_SWITCH;
        double mDensity;
        double mConvectiveVelocityNorm;
        double mEffectiveViscosity;
        double mDeltaTime;
        double mDynamicTau;
        double mTauOne;
        double mTauTwo;
        double mElementSize;
        double mMassProjection;
        double mDynamicViscosity;
        double mPressure;
        double mVelocityDotNabla;

        array_1d<double, 3> mBodyForce;
        array_1d<double, 3> mVelocity;
        array_1d<double, 3> mMeshVelocity;
        array_1d<double, 3> mRelaxedAcceleration;
        array_1d<double, 3> mConvectiveVelocity;
        array_1d<double, 3> mMomentumProjection;
        array_1d<double, 3> mPressureGradient;
        array_1d<double, 3> mEffectiveVelocityDotVelocityGradient;

        BoundedVector<double, TNumNodes> mConvectiveVelocityDotDnDx;
        BoundedVector<double, TNumNodes> mRelaxedAccelerationDotDnDx;
        BoundedVector<double, TNumNodes> mEffectiveVelocityDotVelocityGradientDotShapeGradient;
        BoundedVector<double, TNumNodes> mBodyForceDotDnDx;
        BoundedVector<double, TNumNodes> mMomentumProjectionDotDnDx;
        BoundedVector<double, TNumNodes> mPressureGradientDotDnDx;
        BoundedVector<double, TNumNodes> mNodalPressure;
        BoundedMatrix<double, TNumNodes, TDim> mNodalVelocity;
        BoundedMatrix<double, TNumNodes, TDim> mNodalMeshVelocity;
        BoundedMatrix<double, TNumNodes, TDim> mNodalEffectiveVelocity;
        BoundedMatrix<double, TDim, TDim> mVelocityGradient;
        BoundedMatrix<double, TDim, TDim> mMeshVelocityGradient;
        BoundedMatrix<double, TDim, TDim> mEffectiveVelocityGradient;
        BoundedVector<double, TElementLocalSize> mViscousTermRHSContribution;

        ConstitutiveLaw::Parameters mConstitutiveLawValues;
        BoundedMatrix<double, TStrainSize, TElementLocalSize> mStrainMatrix;
        Vector mStrainRate;
        Vector mShearStress;
        Matrix mC;

        // Sotring this derivatives also in primal data container
        // since all the matrices and vectors needs to be initialized
        // in the heap, therefore to avoid re-reserving memory for each derivative
        ConstitutiveLaw::Parameters mConstitutiveLawValuesDerivative;
        Vector mStrainRateDerivative;
        Vector mShearStressDerivative;
        Matrix mCDerivative;

        ///@}
        ///@name Private Friends
        ///@{

        template<class TDerivativesType, unsigned int TEquationOffset>
        friend class VariableDerivatives;

        ///@}
    };

    ///@}
    ///@name Static Operations
    ///@{

    static array_1d<double, 3> CalculateMomentumProjectionDerivative(
        const IndexType NodeIndex,
        const IndexType DirectionIndex,
        const GeometryType& rGeometry,
        const double Density,
        const array_1d<double, 3>& rBodyForce,
        const array_1d<double, 3>& rPressureGradient,
        const array_1d<double, 3>& rPressureGradientDerivative,
        const array_1d<double, 3>& rEffectiveVelocityDotEffectiveVelocityGradient,
        const array_1d<double, 3>& rEffectiveVelocityDotEffectiveVelocityGradientDerivative)
    {
        // TODO: Implement this for OSS Projection
        return ZeroVector(3);
    }

    static double CalculateMassProjectionDerivative(
        const IndexType NodeIndex,
        const IndexType DirectionIndex,
        const GeometryType& rGeometry)
    {
        // TODO : To be implemented for OSS Projection
        return 0.0;
    }

    static double CalculateNormDerivative(
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

    static void CalculateTau(
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

    static void CalculateTauDerivative(
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
        inv_tau_derivative += -2.0 * c1 * Viscosity * ElementSizeDerivative/ h3;
        inv_tau_derivative += Density * c2 * VelocityNormDerivative / ElementSize;
        inv_tau_derivative += -1.0 * Density * c2 * VelocityNorm * ElementSizeDerivative / h2;
        TauOneDerivative = -1.0 * std::pow(TauOne, 2) * inv_tau_derivative;

        TauTwoDerivative = ViscosityDerivative;
        TauTwoDerivative += c2 * Density * VelocityNormDerivative * ElementSize / c1;
        TauTwoDerivative += c2 * Density * VelocityNorm * ElementSizeDerivative / c1;
    }

    ///@}
};

} // namespace Kratos

#endif // KRATOS_QS_VMS_RESIDUAL_ADJOINT_DERIVATIVES_H
