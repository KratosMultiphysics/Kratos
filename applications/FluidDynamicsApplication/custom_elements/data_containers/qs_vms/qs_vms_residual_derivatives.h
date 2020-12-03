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

// Application includes
#include "custom_constitutive/fluid_constitutive_law.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_elements/data_containers/qs_vms/qs_vms_derivative_utilities.h"

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

    using PropertiesType = typename Element::PropertiesType;

    static constexpr IndexType TBlockSize = TDim + 1;

    /// Size of the strain and stress vectors (in Voigt notation) for the formulation
    constexpr static IndexType TStrainSize = (TDim - 1) * 3; // 3 in 2D, 6 in 3D

    constexpr static IndexType TElementLocalSize = TBlockSize * TNumNodes;

    using ArrayD = array_1d<double, TDim>;

    using VectorN = BoundedVector<double, TNumNodes>;

    using MatrixDD = BoundedMatrix<double, TDim, TDim>;

    using MatrixND = BoundedMatrix<double, TNumNodes, TDim>;

    ///@}
    ///@name Static Operations
    ///@{

    static int Check(
        const GeometryType& rGeometry,
        const ProcessInfo& rProcessInfo);

    static GeometryData::IntegrationMethod GetIntegrationMethod();

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
        ///@name Type Definitions
        ///@{

        constexpr static IndexType TDerivativeDimension = TDerivativesType::TDerivativeDimension;

        ///@}
        ///name@ Life Cycle
        ///@{

        VariableDerivatives(Data& rData) : mrData(rData) {}

        ///@}
        ///@name Operations
        ///@{

        void Initialize(
            Matrix& rOutput,
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

            ArrayD pressure_gradient_derivative;
            BoundedMatrix<double, TDim, TDim> velocity_gradient_derivative;
            FluidCalculationUtilities::EvaluateGradientInPoint(
                r_geometry, rdNdXDerivative,
                std::tie(pressure_gradient_derivative, PRESSURE),
                std::tie(velocity_gradient_derivative, VELOCITY));
            for (IndexType l = 0; l < TDim; ++l) {
                pressure_gradient_derivative[l] += rdNdX(NodeIndex, l) * p_factor;
            }
            row(velocity_gradient_derivative, DirectionIndex) += row(rdNdX, NodeIndex) * u_factor;

            double velocity_dot_nabla_derivative = 0.0;
            for (IndexType a = 0; a < TNumNodes; ++a) {
                for (IndexType i = 0; i < TDim; ++i) {
                    velocity_dot_nabla_derivative += mrData.mNodalVelocity(a, i) * rdNdXDerivative(a, i);
                }
            }

            const ArrayD effective_velocity_derivative_dot_velocity_gradient = prod(mrData.mVelocityGradient, velocity_derivative);
            const ArrayD effective_velocity_dot_velocity_gradient_derivative = prod(velocity_gradient_derivative, mrData.mConvectiveVelocity);
            const VectorN convective_velocity_derivative_dot_dn_dx = prod(rdNdX, velocity_derivative);
            const VectorN relaxed_acceleration_dot_dn_dx_derivative = prod(rdNdXDerivative, mrData.mRelaxedAcceleration);
            const VectorN convective_velocity_dot_dn_dx_derivative = prod(rdNdXDerivative, mrData.mConvectiveVelocity);
            const VectorN effective_velocity_derivative_dot_velocity_gradient_dot_shape_gradient = prod(rdNdX, effective_velocity_derivative_dot_velocity_gradient);
            const VectorN effective_velocity_dot_velocity_gradient_derivative_dot_shape_gradient = prod(rdNdX, effective_velocity_dot_velocity_gradient_derivative);
            const VectorN effective_velocity_dot_velocity_gradient_dot_shape_gradient_derivative = prod(rdNdXDerivative, mrData.mEffectiveVelocityDotVelocityGradient);
            const VectorN pressure_gradient_derivative_dot_shape_gradient = prod(rdNdX, pressure_gradient_derivative);
            const VectorN pressure_gradient_dot_shape_gradient_derivative = prod(rdNdXDerivative, mrData.mPressureGradient);

            // TODO: Needs implementation for OSS projections
            const ArrayD momentum_projection_derivative = ZeroVector(TDim);
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

    class SecondDerivatives
    {
    public:
        ///@name Type Definitions
        ///@{

        ///@}
        ///@name Life Cycle
        ///@{

        SecondDerivatives(
            const Element& rElement,
            FluidConstitutiveLaw& rFluidConstitutiveLaw);

        ///@}
        ///@name Operations
        ///@{

        void Initialize(
            Matrix& rOutput,
            const ProcessInfo& rProcessInfo);

        void AddResidualDerivativeContributions(
            Matrix& rOutput,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

        ///@}
    private:
        ///@name Private Members
        ///@{

        const Element& mrElement;
        FluidConstitutiveLaw& mrFluidConstitutiveLaw;

        IndexType mBlockSize;

        double mDensity;
        double mDynamicViscosity;
        double mElementSize;
        double mDynamicTau;
        double mDeltaTime;

        BoundedMatrix<double, TNumNodes, TDim> mNodalVelocity;

        ConstitutiveLaw::Parameters mConstitutiveLawValues;
        BoundedMatrix<double, TStrainSize, TElementLocalSize> mStrainMatrix;
        Vector mStrainRate;
        Vector mShearStress;
        Matrix mC;

        ///@}

    };

    class Data
    {
    public:
        ///@name Life Cycle
        ///@{

        Data(
            const Element& rElement,
            FluidConstitutiveLaw& rFluidConstitutiveLaw);

        ///@}
        ///@name Operations
        ///@{

        void Initialize(const ProcessInfo& rProcessInfo);

        void CalculateGaussPointData(
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

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

        ArrayD mBodyForce;
        ArrayD mVelocity;
        ArrayD mMeshVelocity;
        ArrayD mRelaxedAcceleration;
        ArrayD mConvectiveVelocity;
        ArrayD mMomentumProjection;
        ArrayD mPressureGradient;
        ArrayD mEffectiveVelocityDotVelocityGradient;

        VectorN mConvectiveVelocityDotDnDx;
        VectorN mRelaxedAccelerationDotDnDx;
        VectorN mEffectiveVelocityDotVelocityGradientDotShapeGradient;
        VectorN mBodyForceDotDnDx;
        VectorN mMomentumProjectionDotDnDx;
        VectorN mPressureGradientDotDnDx;
        VectorN mNodalPressure;
        MatrixND mNodalVelocity;
        MatrixND mNodalMeshVelocity;
        MatrixND mNodalEffectiveVelocity;
        MatrixDD mVelocityGradient;
        MatrixDD mMeshVelocityGradient;
        MatrixDD mEffectiveVelocityGradient;
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

    static double CalculateNormDerivative(
        const double ValueNorm,
        const Vector& Value,
        const Vector& ValueDerivative);

    static void CalculateTau(
        double& TauOne,
        double& TauTwo,
        const double ElementSize,
        const double Density,
        const double Viscosity,
        const double VelocityNorm,
        const double DynamicTau,
        const double DeltaTime);

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
        const double VelocityNormDerivative);

    static void InitializeConstitutiveLaw(
        ConstitutiveLaw::Parameters& rParameters,
        Vector& rStrainVector,
        Vector& rStressVector,
        Matrix& rConstitutiveMatrix,
        const GeometryType& rGeometry,
        const PropertiesType& rProperties,
        const ProcessInfo& rProcessInfo);

    ///@}
};

} // namespace Kratos

#endif // KRATOS_QS_VMS_RESIDUAL_ADJOINT_DERIVATIVES_H
