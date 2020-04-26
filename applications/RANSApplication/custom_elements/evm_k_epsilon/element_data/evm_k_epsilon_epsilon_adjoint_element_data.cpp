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
#include "evm_k_epsilon_epsilon_adjoint_element_data.h"

namespace Kratos
{
namespace EvmKEpsilonAdjointElementDataUtilities
{
template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& EpsilonAdjointElementData<TDim, TNumNodes>::GetAdjointScalarVariable()
{
    return RANS_SCALAR_2_ADJOINT_1;
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& EpsilonAdjointElementData<TDim, TNumNodes>::GetAdjointScalarRateVariable()
{
    return RANS_SCALAR_2_ADJOINT_3;
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::CalculateGaussPointData(
    const Vector& rShapeFunctions, const Matrix& rShapeFunctionDerivatives, const int Step)
{
    KRATOS_TRY

    const double tke = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), TURBULENT_KINETIC_ENERGY, rShapeFunctions, Step);

    this->mKinematicViscosity = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), KINEMATIC_VISCOSITY, rShapeFunctions, Step);
    this->mTurbulentKinematicViscosity = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), TURBULENT_VISCOSITY, rShapeFunctions, Step);
    this->mGamma = EvmKEpsilonElementDataUtilities::CalculateGamma(
        this->mCmu, tke, this->mTurbulentKinematicViscosity);
    this->mVelocityDivergence = RansCalculationUtilities::GetDivergence(
        this->GetGeometry(), VELOCITY, rShapeFunctionDerivatives);

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
void EpsilonAdjointElementData<TDim, TNumNodes>::Check(const GeometryType& rGeometry,
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
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_2_ADJOINT_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_2_ADJOINT_3, r_node);

        KRATOS_CHECK_DOF_IN_NODE(RANS_SCALAR_2_ADJOINT_1, r_node);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Variable<double>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY

    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        noalias(rOutput) = this->mGaussTurbulentKinematicViscositySensitivitiesK *
                           this->mInvEpsilonSigma;
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        noalias(rOutput) = this->mGaussTurbulentKinematicViscositySensitivitiesEpsilon *
                           this->mInvEpsilonSigma;
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const Variable<array_1d<double, 3>>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
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
void EpsilonAdjointElementData<TDim, TNumNodes>::CalculateReactionTermDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Variable<double>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY

    const double reaction =
        this->CalculateReactionTerm(rShapeFunctions, rShapeFunctionDerivatives);

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

    noalias(rOutput) = rOutput * this->mC2;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::CalculateReactionTermDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const Variable<array_1d<double, 3>>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY
    const double reaction =
        this->CalculateReactionTerm(rShapeFunctions, rShapeFunctionDerivatives);
    rOutput.clear();

    if (rDerivativeVariable == VELOCITY)
    {
        if (reaction > 0.0)
        {
            noalias(rOutput) = rShapeFunctionDerivatives * (this->mC1 * 2.0 / 3.0);
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
void EpsilonAdjointElementData<TDim, TNumNodes>::CalculateSourceTermDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Variable<double>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY

    const double production_term =
        EvmKEpsilonElementDataUtilities::CalculateSourceTerm<TDim>(
            this->mVelocityGradient, this->mTurbulentKinematicViscosity);
    BoundedVector<double, TNumNodes> production_sensitivities;

    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        EvmKEpsilonAdjointElementDataUtilities::CalculateProductionScalarSensitivities<TNumNodes>(
            production_sensitivities, this->mTurbulentKinematicViscosity,
            production_term, this->mGaussTurbulentKinematicViscositySensitivitiesK);

        EvmKEpsilonAdjointElementDataUtilities::CalculateGaussThetaTkeSensitivity<TNumNodes>(
            rOutput, this->mCmu, this->mGamma, this->mTurbulentKinematicViscosity,
            this->mGaussTurbulentKinematicViscositySensitivitiesK, rShapeFunctions);
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        EvmKEpsilonAdjointElementDataUtilities::CalculateProductionScalarSensitivities<TNumNodes>(
            production_sensitivities, this->mTurbulentKinematicViscosity, production_term,
            this->mGaussTurbulentKinematicViscositySensitivitiesEpsilon);

        EvmKEpsilonAdjointElementDataUtilities::CalculateGaussThetaEpsilonSensitivity<TNumNodes>(
            rOutput, this->mGamma, this->mTurbulentKinematicViscosity,
            this->mGaussTurbulentKinematicViscositySensitivitiesEpsilon);
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name();
    }

    noalias(rOutput) = rOutput * production_term + production_sensitivities * this->mGamma;
    noalias(rOutput) = rOutput * this->mC1;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::CalculateSourceTermDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const Variable<array_1d<double, 3>>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
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

        noalias(rOutput) = rOutput * (this->mC1 * this->mGamma);
    }
    else
    {
        KRATOS_ERROR << "Unsupported partial derivative variable "
                     << rDerivativeVariable.Name();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
double EpsilonAdjointElementData<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityShapeSensitivity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ShapeParameter& rShapeDerivative,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
double EpsilonAdjointElementData<TDim, TNumNodes>::CalculateReactionTermShapeSensitivity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ShapeParameter& rShapeDerivative,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    const double reaction =
        this->CalculateReactionTerm(rShapeFunctions, rShapeFunctionDerivatives);

    if (reaction > 0.0)
    {
        return (this->mC1 * 2.0 / 3.0) * RansCalculationUtilities::GetDivergence(
                                             this->GetGeometry(), VELOCITY, rDN_Dx_deriv);
    }
    else
    {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double EpsilonAdjointElementData<TDim, TNumNodes>::CalculateSourceTermShapeSensitivity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ShapeParameter& rShapeDerivative,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    const double production_term =
        EvmKEpsilonElementDataUtilities::CalculateSourceTerm<TDim>(
            this->mVelocityGradient, this->mTurbulentKinematicViscosity);

    return EvmKEpsilonAdjointElementDataUtilities::CalculateProductionShapeSensitivities<TDim, TNumNodes>(
               this->mTurbulentKinematicViscosity, 0.0, production_term,
               this->mNodalVelocity, rShapeFunctionDerivatives, rDN_Dx_deriv) *
           this->mC1 * this->mGamma;
}

// template instantiations
template class EpsilonAdjointElementData<2, 3>;
template class EpsilonAdjointElementData<3, 4>;

} // namespace EvmKEpsilonAdjointElementDataUtilities
} // namespace Kratos