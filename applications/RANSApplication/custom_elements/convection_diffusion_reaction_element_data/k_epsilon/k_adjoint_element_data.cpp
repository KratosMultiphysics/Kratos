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
#include "custom_elements/convection_diffusion_reaction_element_data/k_epsilon/adjoint_element_data_utilities.h"
#include "custom_elements/convection_diffusion_reaction_element_data/k_epsilon/element_data_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "k_adjoint_element_data.h"

namespace Kratos
{
///@name Kratos Classes
///@{
namespace KEpsilonElementData
{
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
void KAdjointElementData<TDim, TNumNodes>::CalculateGaussPointData(
    const Vector& rShapeFunctions, const Matrix& rShapeFunctionDerivatives, const int Step)
{
    KRATOS_TRY

    PrimalBaseType::CalculateGaussPointData(rShapeFunctions, rShapeFunctionDerivatives, Step);

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

    KEpsilonElementData::CalculateNodalNutTKESensitivity<TNumNodes>(
        this->mGaussTurbulentKinematicViscositySensitivitiesK, this->mCmu,
        nodal_tke, nodal_epsilon);
    KEpsilonElementData::CalculateNodalNutEpsilonSensitivity<TNumNodes>(
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
void KAdjointElementData<TDim, TNumNodes>::CalculateEffectiveVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const Variable<double>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::CalculateEffectiveVelocityDerivatives(
    BoundedMatrix<double, TNumNodes * TDim, TDim>& rOutput,
    const Variable<array_1d<double, 3>>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();

    for (IndexType c = 0; c < TNumNodes; ++c)
    {
        for (IndexType k = 0; k < TDim; ++k)
        {
            rOutput(c * TDim + k, k) = rShapeFunctions[c];
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, 3> KAdjointElementData<TDim, TNumNodes>::CalculateEffectiveVelocityShapeDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    return ZeroVector(3);
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Variable<double>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
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
double KAdjointElementData<TDim, TNumNodes>::CalculateEffectiveKinematicViscosityShapeDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::CalculateReactionTermDerivatives(
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
            KEpsilonElementData::CalculateGaussGammaTkeSensitivity<TNumNodes>(
                rOutput, this->mCmu, this->mGamma, this->mTurbulentKinematicViscosity,
                this->mGaussTurbulentKinematicViscositySensitivitiesK, rShapeFunctions);
        }
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        if (reaction > 0.0)
        {
            KEpsilonElementData::CalculateGaussGammaEpsilonSensitivity<TNumNodes>(
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
double KAdjointElementData<TDim, TNumNodes>::CalculateReactionTermShapeDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    const double reaction =
        this->CalculateReactionTerm(rShapeFunctions, rShapeFunctionDerivatives);

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
void KAdjointElementData<TDim, TNumNodes>::CalculateSourceTermDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Variable<double>& rDerivativeVariable,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY

    const double production_term = KEpsilonElementData::CalculateSourceTerm<TDim>(
        this->mVelocityGradient, this->mTurbulentKinematicViscosity);

    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        KEpsilonElementData::CalculateProductionScalarSensitivities<TNumNodes>(
            rOutput, this->mTurbulentKinematicViscosity, production_term,
            this->mGaussTurbulentKinematicViscositySensitivitiesK);
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        KEpsilonElementData::CalculateProductionScalarSensitivities<TNumNodes>(
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
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY

    const double production_term = KEpsilonElementData::CalculateSourceTerm<TDim>(
        this->mVelocityGradient, this->mTurbulentKinematicViscosity);

    if (rDerivativeVariable == VELOCITY)
    {
        KEpsilonElementData::CalculateProductionVelocitySensitivities<TDim, TNumNodes>(
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
double KAdjointElementData<TDim, TNumNodes>::CalculateSourceTermShapeDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    const double production_term = KEpsilonElementData::CalculateSourceTerm<TDim>(
        this->mVelocityGradient, this->mTurbulentKinematicViscosity);

    return KEpsilonElementData::CalculateProductionShapeSensitivities<TDim, TNumNodes>(
        this->mTurbulentKinematicViscosity, 0.0, production_term,
        this->mNodalVelocity, rShapeFunctionDerivatives, rDN_Dx_deriv);
}

// template instantiations
template class KAdjointElementData<2, 3>;
template class KAdjointElementData<3, 4>;

} // namespace KEpsilonElementData
///@}
} // namespace Kratos
