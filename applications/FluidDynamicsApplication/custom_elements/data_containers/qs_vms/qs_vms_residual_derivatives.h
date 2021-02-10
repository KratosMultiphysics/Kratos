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
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "utilities/time_discretization.h"

// Application includes
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

    using VectorF = BoundedVector<double, TElementLocalSize>;

    using MatrixDD = BoundedMatrix<double, TDim, TDim>;

    using MatrixND = BoundedMatrix<double, TNumNodes, TDim>;

    ///@}
    ///@name Static Operations
    ///@{

    static void Check(
        const Element& rElement,
        const ProcessInfo& rProcessInfo);

    static GeometryData::IntegrationMethod GetIntegrationMethod();

    ///@}
    ///@name Class Declarations
    ///@{

    class Data;

    ///@}
    ///@name Classes
    ///@{

    class ResidualsContributions
    {
    public:
        ///@name Life Cycle
        ///@{

        ResidualsContributions(Data& rData);

        ///@}
        ///@name Operations
        ///@{

        void AddGaussPointResidualsContributions(
            VectorF& rResidual,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

        ///@}

    private:
        ///@name Private Members
        ///@{

        Data& mrData;

        ///@}
        ///@name Private Operations
        ///@{

        void AddViscousTerms(
            VectorF& rResidual,
            const double W) const;

        ///@}
    };

    template<class TDerivativesType>
    class VariableDerivatives
    {
    public:
        ///@name Type Definitions
        ///@{

        constexpr static IndexType TDerivativeDimension = TDerivativesType::TDerivativeDimension;

        ///@}
        ///name@ Life Cycle
        ///@{

        VariableDerivatives(
            Data& rData)
            : mrData(rData),
              mResidualWeightDerivativeContributions(rData)
        {}

        ///@}
        ///@name Operations
        ///@{

        void CalculateGaussPointResidualsDerivativeContributions(
            VectorF& rResidualDerivative,
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

            // compute viscous term derivative
            double effective_viscosity_derivative;

            // calculate derivative contributions w.r.t. current derivative variable. Derivative variables are
            // assumed be independent of each other, so no cross derivative terms are there. Hence
            // it is only sufficient to derrive w.r.t. current derivative variable.
            mrData.mrConstitutiveLaw.CalculateDerivative(mrData.mConstitutiveLawValues, EFFECTIVE_VISCOSITY, derivatives_type.GetDerivativeVariable(), effective_viscosity_derivative);

            // calculate derivative contributions w.r.t. its dependent variable gradients. Dependent variable gradients
            // may be dependent of each other. therefore it is required to calculate derivatives w.r.t. gradients of
            // all dependent variables.
            double effective_viscosity_derivative_value;
            ArrayD derivative_variable_gradient;
            const auto& r_effective_viscosity_dependent_variables = derivatives_type.GetEffectiveViscosityDependentVariables();
            for (const auto& r_effective_viscosity_dependent_variable : r_effective_viscosity_dependent_variables) {
                // these variables always needs to be scalars (eg. VELOCITY_X)
                const auto& r_derivative_variable = r_effective_viscosity_dependent_variable.GetVariable();
                FluidCalculationUtilities::EvaluateGradientInPoint(mrData.mrElement.GetGeometry(), rdNdXDerivative, std::tie(derivative_variable_gradient, r_derivative_variable));

                // this is a list of gradient component variables. This list also should only contain scalar variables (eg. VELOCITY_GRADIENT_TENSOR_XX)
                const auto& r_effective_viscosity_dependent_variable_gradient_component_list = r_effective_viscosity_dependent_variable.GetVariableGradientComponents();
                for (IndexType i = 0; i < TDim; ++i) {
                    mrData.mrConstitutiveLaw.CalculateDerivative(mrData.mConstitutiveLawValues, EFFECTIVE_VISCOSITY, *r_effective_viscosity_dependent_variable_gradient_component_list[i], effective_viscosity_derivative_value);
                    effective_viscosity_derivative += effective_viscosity_derivative_value * (rdNdX(NodeIndex, i) * (r_derivative_variable ==  derivatives_type.GetDerivativeVariable()));
                    effective_viscosity_derivative += effective_viscosity_derivative_value * (derivative_variable_gradient[i]);
                }
            }

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
                const IndexType row = a * TBlockSize;

                double forcing_derivative = 0.0;
                for (IndexType i = 0; i < TDim; ++i) {

                    double value = 0.0;

                    // Adding RHS derivative terms
                    value += mrData.mDensity * W * tau_one_derivative * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mBodyForce[i];
                    value += mrData.mDensity * W * mrData.mTauOne * convective_velocity_derivative_dot_dn_dx[a] * mrData.mBodyForce[i];
                    value += mrData.mDensity * W * mrData.mTauOne * convective_velocity_dot_dn_dx_derivative[a] * mrData.mBodyForce[i];

                    value -= mrData.mDensity * W * tau_one_derivative * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mMomentumProjection[i];
                    value -= mrData.mDensity * W * mrData.mTauOne * convective_velocity_derivative_dot_dn_dx[a] * mrData.mMomentumProjection[i];
                    value -= mrData.mDensity * W * mrData.mTauOne * convective_velocity_dot_dn_dx_derivative[a] * mrData.mMomentumProjection[i];
                    value -= mrData.mDensity * W * mrData.mTauOne * mrData.mConvectiveVelocityDotDnDx[a] * momentum_projection_derivative[i];

                    value -= W * tau_two_derivative * rdNdX(a, i) * mrData.mMassProjection;
                    value -= W * mrData.mTauTwo * rdNdXDerivative(a, i) * mrData.mMassProjection;
                    value -= W * mrData.mTauTwo * rdNdX(a, i) * mass_projection_derivative;

                    forcing_derivative += rdNdXDerivative(a, i) * mrData.mBodyForce[i];
                    forcing_derivative -= rdNdXDerivative(a, i) * mrData.mMomentumProjection[i];
                    forcing_derivative -= rdNdX(a, i) * momentum_projection_derivative[i];

                    // Adding LHS derivative terms
                    value -= W * mrData.mDensity * rN[a] * effective_velocity_derivative_dot_velocity_gradient[i];
                    value -= W * mrData.mDensity * rN[a] * effective_velocity_dot_velocity_gradient_derivative[i];

                    value -= W * mrData.mDensity * convective_velocity_derivative_dot_dn_dx[a] * mrData.mTauOne * mrData.mDensity * mrData.mEffectiveVelocityDotVelocityGradient[i];
                    value -= W * mrData.mDensity * convective_velocity_dot_dn_dx_derivative[a] * mrData.mTauOne * mrData.mDensity * mrData.mEffectiveVelocityDotVelocityGradient[i];
                    value -= W * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * tau_one_derivative * mrData.mDensity * mrData.mEffectiveVelocityDotVelocityGradient[i];
                    value -= W * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mTauOne * mrData.mDensity * effective_velocity_derivative_dot_velocity_gradient[i];
                    value -= W * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mTauOne * mrData.mDensity * effective_velocity_dot_velocity_gradient_derivative[i];

                    value -= W * tau_one_derivative * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mPressureGradient[i];
                    value -= W * mrData.mTauOne * mrData.mDensity * convective_velocity_derivative_dot_dn_dx[a] * mrData.mPressureGradient[i];
                    value -= W * mrData.mTauOne * mrData.mDensity * convective_velocity_dot_dn_dx_derivative[a] * mrData.mPressureGradient[i];
                    value -= W * mrData.mTauOne * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * pressure_gradient_derivative[i];

                    value += W * rdNdXDerivative(a, i) * mrData.mPressure;
                    value += W * rdNdX(a, i) * rN[NodeIndex] * p_factor;

                    value -= W * tau_two_derivative * rdNdX(a, i) * mrData.mVelocityDotNabla;
                    value -= W * mrData.mTauTwo * rdNdXDerivative(a, i) * mrData.mVelocityDotNabla;
                    value -= W * mrData.mTauTwo * rdNdX(a, i) * velocity_dot_nabla_derivative;
                    value -= W * mrData.mTauTwo * rdNdX(a, i) * rdNdX(NodeIndex, DirectionIndex) * u_factor;

                    // Adding Mass term derivatives

                    value -= W * tau_one_derivative * mrData.mDensity * mrData.mDensity * mrData.mConvectiveVelocityDotDnDx[a] * mrData.mRelaxedAcceleration[i];
                    value -= W * mrData.mTauOne * mrData.mDensity * mrData.mDensity * convective_velocity_derivative_dot_dn_dx[a] * mrData.mRelaxedAcceleration[i];
                    value -= W * mrData.mTauOne * mrData.mDensity * mrData.mDensity * convective_velocity_dot_dn_dx_derivative[a] * mrData.mRelaxedAcceleration[i];

                    rResidualDerivative[row + i] += value;
                }

                double value = 0.0;

                const double forcing = mrData.mBodyForceDotDnDx[a] - mrData.mMomentumProjectionDotDnDx[a];
                value += W * tau_one_derivative * forcing;
                value += W * mrData.mTauOne * forcing_derivative;

                value -= W * tau_one_derivative * mrData.mDensity * mrData.mEffectiveVelocityDotVelocityGradientDotShapeGradient[a];
                value -= W * mrData.mTauOne * mrData.mDensity * effective_velocity_derivative_dot_velocity_gradient_dot_shape_gradient[a];
                value -= W * mrData.mTauOne * mrData.mDensity * effective_velocity_dot_velocity_gradient_derivative_dot_shape_gradient[a];
                value -= W * mrData.mTauOne * mrData.mDensity * effective_velocity_dot_velocity_gradient_dot_shape_gradient_derivative[a];

                value -= W * rN[a] * velocity_dot_nabla_derivative;
                value -= W * rN[a] * rdNdX(NodeIndex, DirectionIndex) * u_factor;

                value -= W * tau_one_derivative * mrData.mPressureGradientDotDnDx[a];
                value -= W * mrData.mTauOne * pressure_gradient_derivative_dot_shape_gradient[a];
                value -= W * mrData.mTauOne * pressure_gradient_dot_shape_gradient_derivative[a];

                // Adding mass term derivatives
                value -=  W * tau_one_derivative * mrData.mDensity * mrData.mRelaxedAccelerationDotDnDx[a];
                value -=  W * mrData.mTauOne * mrData.mDensity * relaxed_acceleration_dot_dn_dx_derivative[a];

                rResidualDerivative[row + TDim] += value;
            }

            mResidualWeightDerivativeContributions.AddGaussPointResidualsContributions(
                rResidualDerivative, WDerivative, rN, rdNdX);

            // calculate shear stress derivative.
            const auto& r_strain_rate_variables = QSVMSDerivativeUtilities<TDim>::GetStrainRateVariables();
            Vector value;
            mrData.mShearStressDerivative.clear();

            // calculate shear stress derivatives w.r.t. strain rate
            for (IndexType i = 0; i < TStrainSize; ++i) {
                mrData.mrConstitutiveLaw.CalculateDerivative(mrData.mConstitutiveLawValues, CAUCHY_STRESS_VECTOR, *r_strain_rate_variables[i], value);
                noalias(mrData.mShearStressDerivative) += value * mrData.mStrainRateDerivative[i];
            }
            // calculate shear stress derivatives w.r.t. effective viscosity
            mrData.mrConstitutiveLaw.CalculateDerivative(mrData.mConstitutiveLawValues, CAUCHY_STRESS_VECTOR, EFFECTIVE_VISCOSITY, value);
            noalias(mrData.mShearStressDerivative) += value * effective_viscosity_derivative;

            this->AddViscousDerivative(rResidualDerivative, NodeIndex,
                                       DirectionIndex, W, rN, rdNdX, WDerivative,
                                       DetJDerivative, rdNdXDerivative);
        }

        ///@}

    private:
        ///@name Private Members
        ///@{

        Data& mrData;
        ResidualsContributions mResidualWeightDerivativeContributions;

        ///@}
        ///@name Private Operations
        ///@{

        void AddViscousDerivative(
            VectorF& rResidualDerivative,
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

            const VectorF& rhs_contribution_derivative =
                prod(trans(mrData.mStrainMatrix), mrData.mShearStressDerivative) +
                prod(trans(strain_matrix_derivative), mrData.mShearStress);

            noalias(rResidualDerivative) -= rhs_contribution_derivative * W;
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
            Data& rData);

        ///@}
        ///@name Operations
        ///@{

        void CalculateGaussPointResidualsDerivativeContributions(
            VectorF& rResidualDerivative,
            const int NodeIndex,
            const int DirectionIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

        ///@}
    private:
        ///@name Private Members
        ///@{

        Data& mrData;

        ///@}

    };

    class Data
    {
    public:
        ///@name Life Cycle
        ///@{

        Data(
            const Element& rElement,
            ConstitutiveLaw& rConstitutiveLaw,
            const ProcessInfo& rProcessInfo);

        ///@}
        ///@name Operations
        ///@{

        void CalculateGaussPointData(
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

        ///@}
    private:
        ///@name Private Members
        ///@{

        const Element& mrElement;
        ConstitutiveLaw& mrConstitutiveLaw;

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
        Vector mStrainRateDerivative;
        Vector mShearStressDerivative;

        ///@}
        ///@name Private Friends
        ///@{

        template<class TDerivativesType>
        friend class VariableDerivatives;
        friend class SecondDerivatives;
        friend class ResidualsContributions;

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
