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

#pragma once

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

/**
 * @brief Computes the QSVMS residual derivatives
 *
 * This class holds sub-classes which are required to compute QSVMS residual derivatives
 * with respect to given variables.
 *
 * @tparam TDim
 * @tparam TNumNodes
 */
template <unsigned int TDim, unsigned int TNumNodes>
class QSVMSResidualDerivatives
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using PropertiesType = typename Element::PropertiesType;

    static constexpr IndexType TBlockSize = TDim + 1;

    /// Size of the strain and stress vectors (in Voigt notation) for the formulation
    constexpr static IndexType TStrainSize = (TDim - 1) * 3; // 3 in 2D, 6 in 3D

    constexpr static IndexType TElementLocalSize = TBlockSize * TNumNodes;

    constexpr static IndexType TNN = TNumNodes;

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

    /**
     * @brief The data container for the QSVMS residuals
     *
     * This data container is used for all the derivative computations. So, this
     * is used to store common data which are required for derivative computations
     * between different variables.
     *
     */
    class QSVMSResidualData;

    ///@}
    ///@name Classes
    ///@{

    /**
     * @brief Computes residual contributions
     *
     * This class is used to compute the gauss point residuals.
     *
     */
    class ResidualsContributions
    {
    public:
        ///@name Type definitions
        ///@{

        constexpr static IndexType NumNodes = TNumNodes;

        constexpr static IndexType BlockSize = TBlockSize;

        ///@}
        ///@name Operations
        ///@{

        void AddGaussPointResidualsContributions(
            VectorF& rResidual,
            QSVMSResidualData& rData,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX) const;

        ///@}

    private:
        ///@name Private Operations
        ///@{

        void static AddViscousTerms(
            QSVMSResidualData& rData,
            VectorF& rResidual,
            const double W);

        ///@}
    };

    /**
     * @brief Computes QS VMS residual derivative residuals for given variable
     *
     * This class is capable of computing QS VMS redisual derivative with respect to
     * given derivative variable (information is provided via TDerivativesType)
     *
     * @tparam TDerivativesType     Derivative computation class which holds information on how to compute coefficient derivatives
     */
    template<class TDerivativesType>
    class VariableDerivatives
    {
    public:
        ///@name Type Definitions
        ///@{

        constexpr static IndexType NumNodes = TNumNodes;

        constexpr static IndexType BlockSize = TBlockSize;

        constexpr static IndexType TDerivativeDimension = TDerivativesType::TDerivativeDimension;

        constexpr static IndexType ComponentIndex = TDerivativesType::ComponentIndex;

        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief Computes gauss point residual derivatives for given NodeIndex
         *
         * This method computes gauss point QS VMS residual derivatives w.r.t. the variable information in TDerivativesType
         * for the NodeIndex.
         *
         * @param rResidualDerivative           Output derivatives vector of the Residual
         * @param rData                         Data holder, which holds all the intermediate values computed from the primal solution which are common to all the derivatives
         * @param NodeIndex                     Index of the node in the geometry for which the derivatives are computed for.
         * @param W                             Gauss point weight
         * @param rN                            Gauss point shape function values vector
         * @param rdNdX                         Gauss point shape gradient matrix. (rows -> nodes, columns -> coordinates index)
         * @param WDerivative                   Gauss point weight derivative w.r.t. chosen variable from TDerivativesType
         * @param DetJDerivative                Jacobian determinant derivative w.r.t. chosen variable from TDerivativesType
         * @param rdNdXDerivative               Gauss point shape gradient's derivatives w.r.t. chosen variable from TDerivativesType
         * @param MassTermsDerivativesWeight    This decides whether Mass term derivatives needs to be computed or not (either 1.0 or 0.0)
         */
        void CalculateGaussPointResidualsDerivativeContributions(
            VectorF& rResidualDerivative,
            QSVMSResidualData& rData,
            const int NodeIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative,
            const double MassTermsDerivativesWeight = 1.0) const
        {
            rResidualDerivative.clear();

            constexpr double u_factor = TDerivativesType::VelocityDerivativeFactor;
            constexpr double p_factor = TDerivativesType::PressureDerivativeFactor;

            const auto& r_geometry = rData.mpElement->GetGeometry();

            const TDerivativesType derivatives_type(
                NodeIndex, r_geometry, W, rN, rdNdX,
                WDerivative, DetJDerivative, rdNdXDerivative);

            const auto& velocity_derivative = derivatives_type.CalculateEffectiveVelocityDerivative(rData.mVelocity);
            const auto& element_length_derivative = derivatives_type.CalculateElementLengthDerivative(rData.mElementSize);
            derivatives_type.CalculateStrainRateDerivative(rData.mStrainRateDerivative, rData.mNodalVelocity);

            // compute viscous term derivative
            double effective_viscosity_derivative;

            // calculate derivative contributions w.r.t. current derivative variable. Derivative variables are
            // assumed be independent of each other, so no cross derivative terms are there. Hence
            // it is only sufficient to derive w.r.t. current derivative variable.
            rData.mpConstitutiveLaw->CalculateDerivative(rData.mConstitutiveLawValues, EFFECTIVE_VISCOSITY, derivatives_type.GetDerivativeVariable(), effective_viscosity_derivative);
            effective_viscosity_derivative *= rN[NodeIndex];

            // calculate derivative contributions w.r.t. its dependent variable gradients. Dependent variable gradients
            // may be dependent of each other. therefore it is required to calculate derivatives w.r.t. gradients of
            // all dependent variables.
            double effective_viscosity_derivative_value;
            ArrayD derivative_variable_gradient;
            const auto& r_effective_viscosity_dependent_variables = derivatives_type.GetEffectiveViscosityDependentVariables();
            for (const auto& r_effective_viscosity_dependent_variable : r_effective_viscosity_dependent_variables) {
                // these variables always needs to be scalars (eg. VELOCITY_X)
                const auto& r_derivative_variable = std::get<0>(r_effective_viscosity_dependent_variable);
                FluidCalculationUtilities::EvaluateGradientInPoint(rData.mpElement->GetGeometry(), rdNdXDerivative, std::tie(derivative_variable_gradient, r_derivative_variable));

                // this is a list of gradient component variables. This list also should only contain scalar variables (eg. VELOCITY_GRADIENT_TENSOR_XX)
                const auto& r_effective_viscosity_dependent_variable_gradient_component_list = std::get<1>(r_effective_viscosity_dependent_variable);
                for (IndexType i = 0; i < TDim; ++i) {
                    rData.mpConstitutiveLaw->CalculateDerivative(rData.mConstitutiveLawValues, EFFECTIVE_VISCOSITY, *r_effective_viscosity_dependent_variable_gradient_component_list[i], effective_viscosity_derivative_value);
                    effective_viscosity_derivative += effective_viscosity_derivative_value * (rdNdX(NodeIndex, i) * (r_derivative_variable ==  derivatives_type.GetDerivativeVariable()));
                    effective_viscosity_derivative += effective_viscosity_derivative_value * (derivative_variable_gradient[i]);
                }
            }

            const double velocity_norm_derivative = CalculateNormDerivative(
                rData.mConvectiveVelocityNorm, rData.mConvectiveVelocity, velocity_derivative);

            double tau_one_derivative, tau_two_derivative;
            CalculateTauDerivative(
                tau_one_derivative, tau_two_derivative, rData.mTauOne,
                rData.mDensity, rData.mDynamicTau, rData.mDeltaTime,
                rData.mElementSize, element_length_derivative,
                rData.mEffectiveViscosity, effective_viscosity_derivative,
                rData.mConvectiveVelocityNorm, velocity_norm_derivative);

            ArrayD pressure_gradient_derivative;
            BoundedMatrix<double, TDim, TDim> velocity_gradient_derivative;
            FluidCalculationUtilities::EvaluateGradientInPoint(
                r_geometry, rdNdXDerivative,
                std::tie(pressure_gradient_derivative, PRESSURE),
                std::tie(velocity_gradient_derivative, VELOCITY));
            for (IndexType l = 0; l < TDim; ++l) {
                pressure_gradient_derivative[l] += rdNdX(NodeIndex, l) * p_factor;
            }
            row(velocity_gradient_derivative, ComponentIndex) += row(rdNdX, NodeIndex) * u_factor;

            double velocity_dot_nabla_derivative = 0.0;
            for (IndexType a = 0; a < TNumNodes; ++a) {
                for (IndexType i = 0; i < TDim; ++i) {
                    velocity_dot_nabla_derivative += rData.mNodalVelocity(a, i) * rdNdXDerivative(a, i);
                }
            }

            const ArrayD effective_velocity_derivative_dot_velocity_gradient = prod(rData.mVelocityGradient, velocity_derivative);
            const ArrayD effective_velocity_dot_velocity_gradient_derivative = prod(velocity_gradient_derivative, rData.mConvectiveVelocity);
            const VectorN convective_velocity_derivative_dot_dn_dx = prod(rdNdX, velocity_derivative);
            const VectorN relaxed_acceleration_dot_dn_dx_derivative = prod(rdNdXDerivative, rData.mRelaxedAcceleration);
            const VectorN convective_velocity_dot_dn_dx_derivative = prod(rdNdXDerivative, rData.mConvectiveVelocity);
            const VectorN effective_velocity_derivative_dot_velocity_gradient_dot_shape_gradient = prod(rdNdX, effective_velocity_derivative_dot_velocity_gradient);
            const VectorN effective_velocity_dot_velocity_gradient_derivative_dot_shape_gradient = prod(rdNdX, effective_velocity_dot_velocity_gradient_derivative);
            const VectorN effective_velocity_dot_velocity_gradient_dot_shape_gradient_derivative = prod(rdNdXDerivative, rData.mEffectiveVelocityDotVelocityGradient);
            const VectorN pressure_gradient_derivative_dot_shape_gradient = prod(rdNdX, pressure_gradient_derivative);
            const VectorN pressure_gradient_dot_shape_gradient_derivative = prod(rdNdXDerivative, rData.mPressureGradient);

            // TODO: Needs implementation for OSS projections
            const ArrayD momentum_projection_derivative = ZeroVector(TDim);
            const double mass_projection_derivative = 0.0;

            for (IndexType a = 0; a < TNumNodes; ++a) {
                const IndexType row = a * TBlockSize;

                double forcing_derivative = 0.0;
                for (IndexType i = 0; i < TDim; ++i) {

                    double value = 0.0;

                    // Adding RHS derivative terms
                    value += rData.mDensity * W * tau_one_derivative * rData.mConvectiveVelocityDotDnDx[a] * rData.mBodyForce[i];
                    value += rData.mDensity * W * rData.mTauOne * convective_velocity_derivative_dot_dn_dx[a] * rData.mBodyForce[i];
                    value += rData.mDensity * W * rData.mTauOne * convective_velocity_dot_dn_dx_derivative[a] * rData.mBodyForce[i];

                    value -= rData.mDensity * W * tau_one_derivative * rData.mConvectiveVelocityDotDnDx[a] * rData.mMomentumProjection[i];
                    value -= rData.mDensity * W * rData.mTauOne * convective_velocity_derivative_dot_dn_dx[a] * rData.mMomentumProjection[i];
                    value -= rData.mDensity * W * rData.mTauOne * convective_velocity_dot_dn_dx_derivative[a] * rData.mMomentumProjection[i];
                    value -= rData.mDensity * W * rData.mTauOne * rData.mConvectiveVelocityDotDnDx[a] * momentum_projection_derivative[i];

                    value -= W * tau_two_derivative * rdNdX(a, i) * rData.mMassProjection;
                    value -= W * rData.mTauTwo * rdNdXDerivative(a, i) * rData.mMassProjection;
                    value -= W * rData.mTauTwo * rdNdX(a, i) * mass_projection_derivative;

                    forcing_derivative += rdNdXDerivative(a, i) * rData.mBodyForce[i];
                    forcing_derivative -= rdNdXDerivative(a, i) * rData.mMomentumProjection[i];
                    forcing_derivative -= rdNdX(a, i) * momentum_projection_derivative[i];

                    // Adding LHS derivative terms
                    value -= W * rData.mDensity * rN[a] * effective_velocity_derivative_dot_velocity_gradient[i];
                    value -= W * rData.mDensity * rN[a] * effective_velocity_dot_velocity_gradient_derivative[i];

                    value -= W * rData.mDensity * convective_velocity_derivative_dot_dn_dx[a] * rData.mTauOne * rData.mDensity * rData.mEffectiveVelocityDotVelocityGradient[i];
                    value -= W * rData.mDensity * convective_velocity_dot_dn_dx_derivative[a] * rData.mTauOne * rData.mDensity * rData.mEffectiveVelocityDotVelocityGradient[i];
                    value -= W * rData.mDensity * rData.mConvectiveVelocityDotDnDx[a] * tau_one_derivative * rData.mDensity * rData.mEffectiveVelocityDotVelocityGradient[i];
                    value -= W * rData.mDensity * rData.mConvectiveVelocityDotDnDx[a] * rData.mTauOne * rData.mDensity * effective_velocity_derivative_dot_velocity_gradient[i];
                    value -= W * rData.mDensity * rData.mConvectiveVelocityDotDnDx[a] * rData.mTauOne * rData.mDensity * effective_velocity_dot_velocity_gradient_derivative[i];

                    value -= W * tau_one_derivative * rData.mDensity * rData.mConvectiveVelocityDotDnDx[a] * rData.mPressureGradient[i];
                    value -= W * rData.mTauOne * rData.mDensity * convective_velocity_derivative_dot_dn_dx[a] * rData.mPressureGradient[i];
                    value -= W * rData.mTauOne * rData.mDensity * convective_velocity_dot_dn_dx_derivative[a] * rData.mPressureGradient[i];
                    value -= W * rData.mTauOne * rData.mDensity * rData.mConvectiveVelocityDotDnDx[a] * pressure_gradient_derivative[i];

                    value += W * rdNdXDerivative(a, i) * rData.mPressure;
                    value += W * rdNdX(a, i) * rN[NodeIndex] * p_factor;

                    value -= W * tau_two_derivative * rdNdX(a, i) * rData.mVelocityDotNabla;
                    value -= W * rData.mTauTwo * rdNdXDerivative(a, i) * rData.mVelocityDotNabla;
                    value -= W * rData.mTauTwo * rdNdX(a, i) * velocity_dot_nabla_derivative;
                    value -= W * rData.mTauTwo * rdNdX(a, i) * rdNdX(NodeIndex, ComponentIndex) * u_factor;

                    // Adding Mass term derivatives

                    value -= W * tau_one_derivative * rData.mDensity * rData.mDensity * rData.mConvectiveVelocityDotDnDx[a] * rData.mRelaxedAcceleration[i] * MassTermsDerivativesWeight;
                    value -= W * rData.mTauOne * rData.mDensity * rData.mDensity * convective_velocity_derivative_dot_dn_dx[a] * rData.mRelaxedAcceleration[i] * MassTermsDerivativesWeight;
                    value -= W * rData.mTauOne * rData.mDensity * rData.mDensity * convective_velocity_dot_dn_dx_derivative[a] * rData.mRelaxedAcceleration[i] * MassTermsDerivativesWeight;

                    rResidualDerivative[row + i] += value;
                }

                double value = 0.0;

                const double forcing = rData.mBodyForceDotDnDx[a] - rData.mMomentumProjectionDotDnDx[a];
                value += W * tau_one_derivative * forcing;
                value += W * rData.mTauOne * forcing_derivative;

                value -= W * tau_one_derivative * rData.mDensity * rData.mEffectiveVelocityDotVelocityGradientDotShapeGradient[a];
                value -= W * rData.mTauOne * rData.mDensity * effective_velocity_derivative_dot_velocity_gradient_dot_shape_gradient[a];
                value -= W * rData.mTauOne * rData.mDensity * effective_velocity_dot_velocity_gradient_derivative_dot_shape_gradient[a];
                value -= W * rData.mTauOne * rData.mDensity * effective_velocity_dot_velocity_gradient_dot_shape_gradient_derivative[a];

                value -= W * rN[a] * velocity_dot_nabla_derivative;
                value -= W * rN[a] * rdNdX(NodeIndex, ComponentIndex) * u_factor;

                value -= W * tau_one_derivative * rData.mPressureGradientDotDnDx[a];
                value -= W * rData.mTauOne * pressure_gradient_derivative_dot_shape_gradient[a];
                value -= W * rData.mTauOne * pressure_gradient_dot_shape_gradient_derivative[a];

                // Adding mass term derivatives
                value -=  W * tau_one_derivative * rData.mDensity * rData.mRelaxedAccelerationDotDnDx[a] * MassTermsDerivativesWeight;
                value -=  W * rData.mTauOne * rData.mDensity * relaxed_acceleration_dot_dn_dx_derivative[a] * MassTermsDerivativesWeight;

                rResidualDerivative[row + TDim] += value;
            }

            mResidualsContributions.AddGaussPointResidualsContributions(
                rResidualDerivative, rData, WDerivative, rN, rdNdX);

            // calculate shear stress derivative.
            const auto& r_strain_rate_variables = QSVMSDerivativeUtilities<TDim>::GetStrainRateVariables();
            Vector value;
            rData.mShearStressDerivative.clear();

            // calculate shear stress derivatives w.r.t. strain rate
            for (IndexType i = 0; i < TStrainSize; ++i) {
                rData.mpConstitutiveLaw->CalculateDerivative(rData.mConstitutiveLawValues, CAUCHY_STRESS_VECTOR, *r_strain_rate_variables[i], value);
                noalias(rData.mShearStressDerivative) += value * rData.mStrainRateDerivative[i];
            }
            // calculate shear stress derivatives w.r.t. effective viscosity
            rData.mpConstitutiveLaw->CalculateDerivative(rData.mConstitutiveLawValues, CAUCHY_STRESS_VECTOR, EFFECTIVE_VISCOSITY, value);
            noalias(rData.mShearStressDerivative) += value * effective_viscosity_derivative;

            AddViscousDerivative(rData, rResidualDerivative, NodeIndex,
                                       W, rN, rdNdX, WDerivative,
                                       DetJDerivative, rdNdXDerivative);
        }

        ///@}

    private:
        ///@name Private members
        ///@{

        ResidualsContributions mResidualsContributions;

        ///@}
        ///@name Private Operations
        ///@{

        /**
         * @brief Adds the viscous term derivative of the residual
         *
         * @param rData                         Data holder, which holds all the intermediate values computed from the primal solution which are common to all the derivatives
         * @param rResidualDerivative           Output derivatives vector of the Residual
         * @param NodeIndex                     Index of the node in the geometry for which the derivatives are computed for.
         * @param W                             Gauss point weight
         * @param rN                            Gauss point shape function values vector
         * @param rdNdX                         Gauss point shape gradient matrix. (rows -> nodes, columns -> coordinates index)
         * @param WDerivative                   Gauss point weight derivative w.r.t. chosen variable from TDerivativesType
         * @param DetJDerivative                Jacobian determinant derivative w.r.t. chosen variable from TDerivativesType
         * @param rdNdXDerivative               auss point shape gradient's derivatives w.r.t. chosen variable from TDerivativesType
         */
        void static AddViscousDerivative(
            QSVMSResidualData& rData,
            VectorF& rResidualDerivative,
            const int NodeIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX,
            const double WDerivative,
            const double DetJDerivative,
            const Matrix& rdNdXDerivative)
        {
            BoundedMatrix<double, TStrainSize, TElementLocalSize> strain_matrix_derivative;
            FluidElementUtilities<TNumNodes>::GetStrainMatrix(rdNdXDerivative, strain_matrix_derivative);

            const VectorF& rhs_contribution_derivative =
                prod(trans(rData.mStrainMatrix), rData.mShearStressDerivative) +
                prod(trans(strain_matrix_derivative), rData.mShearStress);

            noalias(rResidualDerivative) -= rhs_contribution_derivative * W;
        }

        ///@}
    };

    /**
     * @brief Computes second derivatives of the QS VMS residual
     *
     * This class computes second derivatives with respect to given state variables
     *
     * @tparam TComponentIndex
     */
    template<unsigned int TComponentIndex>
    class SecondDerivatives
    {
    public:
        ///@name Type definitions
        ///@{

        constexpr static IndexType NumNodes = TNumNodes;

        constexpr static IndexType BlockSize = TBlockSize;

        ///@}
        ///@name Operations
        ///@{

        void CalculateGaussPointResidualsDerivativeContributions(
            VectorF& rResidualDerivative,
            QSVMSResidualData& rData,
            const int NodeIndex,
            const double W,
            const Vector& rN,
            const Matrix& rdNdX) const;

        ///@}
    };

    class QSVMSResidualData
    {
    public:
        ///@name Type definitions
        ///@{

        static constexpr IndexType TBlockSize = TDim + 1;

        static constexpr IndexType TLNumNodes = TNN;

        static constexpr IndexType TResidualSize = TBlockSize * TLNumNodes;

        ///@}
        ///@name Operations
        ///@{

        void Initialize(
            const Element& rElement,
            ConstitutiveLaw& rConstitutiveLaw,
            const ProcessInfo& rProcessInfo);

        void CalculateGaussPointData(
            const double W,
            const Vector& rN,
            const Matrix& rdNdX);

        ///@}
    private:
        ///@name Private Members
        ///@{

        const Element* mpElement;
        ConstitutiveLaw* mpConstitutiveLaw;

        // Primal data
        int mOSSSwitch;
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

        // Sorting this derivatives also in primal data container
        // since all the matrices and vectors needs to be initialized
        // in the heap, therefore to avoid re-reserving memory for each derivative
        Vector mStrainRateDerivative;
        Vector mShearStressDerivative;

        ///@}
        ///@name Private Friends
        ///@{

        template<class TDerivativesType>
        friend class VariableDerivatives;

        template<unsigned int TComponentIndex>
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