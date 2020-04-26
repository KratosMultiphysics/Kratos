//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/variables.h"

// Application includes
#include "custom_elements/evm_k_epsilon/element_data/evm_k_epsilon_adjoint_element_data_utilities.h"
#include "custom_elements/evm_k_epsilon/element_data/evm_k_epsilon_element_data_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "evm_k_epsilon_k_adjoint_element_data.h"

namespace Kratos
{
namespace EvmKEpsilonAdjointElementDataUtilities
{
template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& KAdjointElementData<TDim, TNumNodes>::GetAdjointScalarVariable()
{
    return RANS_SCALAR_1_ADJOINT_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& KAdjointElementData<TDim, TNumNodes>::GetAdjointScalarRateVariable()
{
    return RANS_SCALAR_1_ADJOINT_3;
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::CalculateGaussPointData(const Vector& rShapeFunctions,
                                                                   const Matrix& rShapeFunctionDerivatives,
                                                                   const ProcessInfo& rCurrentProcessInfo,
                                                                   const int Step)
{
    KRATOS_TRY

    this->mInvTkeSigma = 1.0 / rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];
    this->mCmu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    this->mTurbulentKinematicViscosity = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), TURBULENT_VISCOSITY, rShapeFunctions, Step);
    this->mTurbulentKineticEnergy = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), TURBULENT_KINETIC_ENERGY, rShapeFunctions, Step);
    this->mKinematicViscosity = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), KINEMATIC_VISCOSITY, rShapeFunctions, Step);
    this->mGamma = EvmKEpsilonElementDataUtilities::CalculateGamma(
        this->mCmu, this->mTurbulentKineticEnergy, this->mTurbulentKinematicViscosity);

    this->mVelocityDivergence = RansCalculationUtilities::GetDivergence(
        this->GetGeometry(), VELOCITY, rShapeFunctionDerivatives, Step);

    RansCalculationUtilities::CalculateGradient<TDim>(
        this->mVelocityGradient, this->GetGeometry(), VELOCITY,
        rShapeFunctionDerivatives, Step);

    BoundedVector<double, TNumNodes> nodal_tke;
    BoundedVector<double, TNumNodes> nodal_epsilon;

    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
    {
        const NodeType& r_node = this->GetGeometry()[i_node];
        const array_1d<double, 3>& rVelocity =
            r_node.FastGetSolutionStepValue(VELOCITY, Step);
        for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
        {
            mNodalVelocity(i_node, i_dim) = rVelocity[i_dim];
        }
        nodal_tke[i_node] =
            r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY, Step);
        nodal_epsilon[i_node] =
            r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE, Step);
    }

    EvmKEpsilonAdjointElementDataUtilities::CalculateNodalNutTKESensitivity<TNumNodes>(
        this->mGaussTurbulentKinematicViscositySensitivitiesK, this->mCmu,
        nodal_tke, nodal_epsilon);
    EvmKEpsilonAdjointElementDataUtilities::CalculateNodalNutEpsilonSensitivity<TNumNodes>(
        this->mGaussTurbulentKinematicViscositySensitivitiesEpsilon, this->mCmu,
        nodal_tke, nodal_epsilon);

    RansCalculationUtilities::CalculateGaussSensitivities<TNumNodes>(
        this->mGaussTurbulentKinematicViscositySensitivitiesK,
        this->mGaussTurbulentKinematicViscositySensitivitiesK, rShapeFunctions);
    RansCalculationUtilities::CalculateGaussSensitivities<TNumNodes>(
        this->mGaussTurbulentKinematicViscositySensitivitiesEpsilon,
        this->mGaussTurbulentKinematicViscositySensitivitiesEpsilon, rShapeFunctions);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::Check(const GeometryType& rGeometry,
                                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const int number_of_nodes = rGeometry.PointsNumber();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const NodeType& r_node = rGeometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_1_ADJOINT_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_1_ADJOINT_3, r_node);

        KRATOS_CHECK_DOF_IN_NODE(RANS_SCALAR_1_ADJOINT_1, r_node);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Variable<double>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        noalias(rOutput) =
            this->mGaussTurbulentKinematicViscositySensitivitiesK * this->mInvTkeSigma;
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        noalias(rOutput) =
            this->mGaussTurbulentKinematicViscositySensitivitiesEpsilon * this->mInvTkeSigma;
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const Variable<array_1d<double, 3>>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rDerivativeVariable == VELOCITY)
    {
        rOutput.clear();
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::CalculateReactionTermDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Variable<double>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const double reaction = this->CalculateReactionTerm(
        rShapeFunctions, rShapeFunctionDerivatives, rCurrentProcessInfo);

    rOutput.clear();

    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        if (reaction > 0.0)
        {
            EvmKEpsilonAdjointElementDataUtilities::CalculateGaussThetaTkeSensitivity<TNumNodes>(
                rOutput, this->mCmu, this->mGamma, this->mTurbulentKinematicViscosity,
                this->mGaussTurbulentKinematicViscositySensitivitiesK, rShapeFunctions);
        }
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        if (reaction > 0.0)
        {
            EvmKEpsilonAdjointElementDataUtilities::CalculateGaussThetaEpsilonSensitivity<TNumNodes>(
                rOutput, this->mGamma, this->mTurbulentKinematicViscosity,
                this->mGaussTurbulentKinematicViscositySensitivitiesEpsilon);
        }
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::CalculateReactionTermDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const Variable<array_1d<double, 3>>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    const double reaction = this->CalculateReactionTerm(
        rShapeFunctions, rShapeFunctionDerivatives, rCurrentProcessInfo);
    rOutput.clear();

    if (rDerivativeVariable == VELOCITY)
    {
        if (reaction > 0.0)
        {
            noalias(rOutput) = rShapeFunctionDerivatives * (2.0 / 3.0);
        }
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::CalculateSourceTermDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Variable<double>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const double production_term =
        EvmKEpsilonElementDataUtilities::CalculateSourceTerm<TDim>(
            this->mVelocityGradient, this->mTurbulentKinematicViscosity);

    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        EvmKEpsilonAdjointElementDataUtilities::CalculateProductionScalarSensitivities<TNumNodes>(
            rOutput, this->mTurbulentKinematicViscosity, production_term,
            this->mGaussTurbulentKinematicViscositySensitivitiesK);
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        EvmKEpsilonAdjointElementDataUtilities::CalculateProductionScalarSensitivities<TNumNodes>(
            rOutput, this->mTurbulentKinematicViscosity, production_term,
            this->mGaussTurbulentKinematicViscositySensitivitiesEpsilon);
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::CalculateSourceTermDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const Variable<array_1d<double, 3>>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const double production_term =
        EvmKEpsilonElementDataUtilities::CalculateSourceTerm<TDim>(
            this->mVelocityGradient, this->mTurbulentKinematicViscosity);

    if (rDerivativeVariable == VELOCITY)
    {
        EvmKEpsilonAdjointElementDataUtilities::CalculateProductionVelocitySensitivities<TDim, TNumNodes>(
            rOutput, this->mTurbulentKinematicViscosity, production_term,
            this->mVelocityGradient, rShapeFunctionDerivatives);
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
double KAdjointElementData<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityShapeSensitivity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ShapeParameter& rShapeDerivative,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
    const ProcessInfo& rCurrentProcessInfo) const
{
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
double KAdjointElementData<TDim, TNumNodes>::CalculateReactionTermShapeSensitivity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ShapeParameter& rShapeDerivative,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const double reaction = this->CalculateReactionTerm(
        rShapeFunctions, rShapeFunctionDerivatives, rCurrentProcessInfo);

    if (reaction > 0.0)
    {
        return (2.0 / 3.0) * RansCalculationUtilities::GetDivergence(
                                 this->GetGeometry(), VELOCITY, rDN_Dx_deriv);
    }
    else
    {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double KAdjointElementData<TDim, TNumNodes>::CalculateSourceTermShapeSensitivity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ShapeParameter& rShapeDerivative,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const double production_term =
        EvmKEpsilonElementDataUtilities::CalculateSourceTerm<TDim>(
            this->mVelocityGradient, this->mTurbulentKinematicViscosity);

    return EvmKEpsilonAdjointElementDataUtilities::CalculateProductionShapeSensitivities<TDim, TNumNodes>(
        this->mTurbulentKinematicViscosity, 0.0, production_term,
        this->mNodalVelocity, rShapeFunctionDerivatives, rDN_Dx_deriv);
}

// template instantiations
template class KAdjointElementData<2, 3>;
template class KAdjointElementData<3, 4>;

} // namespace EvmKEpsilonAdjointElementDataUtilities
} // namespace Kratos