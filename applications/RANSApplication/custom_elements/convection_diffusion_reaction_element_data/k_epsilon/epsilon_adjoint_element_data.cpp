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
#include "epsilon_adjoint_element_data.h"

namespace Kratos
{
///@name Kratos Classes
///@{
namespace KEpsilonElementData
{
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
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_2_ADJOINT_1, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_2_ADJOINT_3, r_node);

        KRATOS_CHECK_DOF_IN_NODE(RANS_SCALAR_2_ADJOINT_1, r_node);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::CalculateGaussPointData(
    const Vector& rShapeFunctions, const Matrix& rShapeFunctionDerivatives, const int Step)
{
    KRATOS_TRY

    PrimalBaseType::CalculateGaussPointData(rShapeFunctions, rShapeFunctionDerivatives, Step);

    mReactionTerm = this->CalculateReactionTerm(rShapeFunctions, rShapeFunctionDerivatives);
    mProductionTerm = KEpsilonElementData::CalculateSourceTerm<TDim>(
        this->mVelocityGradient, this->mTurbulentKinematicViscosity);

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
const typename EpsilonAdjointElementData<TDim, TNumNodes>::AdjointBaseType::ScalarDerivative& EpsilonAdjointElementData<TDim, TNumNodes>::GetScalarDerivativeData(
    const Variable<double>& rDerivativeVariable) const
{
    KRATOS_TRY

    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY)
    {
        return mKDerivativeData;
    }
    else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE)
    {
        return mEpsilonDerivativeData;
    }
    else
    {
        KRATOS_ERROR << "Unsupported " << rDerivativeVariable.Name()
                     << " derivative variable for " << this->GetName();
        return mKDerivativeData;
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
const typename EpsilonAdjointElementData<TDim, TNumNodes>::AdjointBaseType::VectorDerivative& EpsilonAdjointElementData<TDim, TNumNodes>::GetVectorDerivativeData(
    const Variable<array_1d<double, 3>>& rDerivativeVariable) const
{
    // Currently only velocity derivatives are required, no checks for
    // variables are done. This is helpfull in creating unit tests
    return mVelocityDerivativeData;
}

template <unsigned int TDim, unsigned int TNumNodes>
const typename EpsilonAdjointElementData<TDim, TNumNodes>::AdjointBaseType::ShapeDerivative& EpsilonAdjointElementData<TDim, TNumNodes>::GetShapeDerivativeData(
    const Variable<array_1d<double, 3>>& rDerivativeVariable) const
{
    // Currently only shape derivatives are required, no checks for
    // variables are done. This is helpfull in creating unit tests
    return mShapeDerivativeData;
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::KDerivativeData::CalculateEffectiveVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::EpsilonDerivativeData::CalculateEffectiveVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::VelocityDerivativeData::CalculateEffectiveVelocityDerivatives(
    BoundedMatrix<double, TNumNodes * TDim, TDim>& rOutput,
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
array_1d<double, 3> EpsilonAdjointElementData<TDim, TNumNodes>::ShapeDerivativeData::CalculateEffectiveVelocityShapeDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    return ZeroVector(3);
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::KDerivativeData::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY

    noalias(rOutput) = mData.mGaussTurbulentKinematicViscositySensitivitiesK *
                       mData.mInvEpsilonSigma;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::EpsilonDerivativeData::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY

    noalias(rOutput) = mData.mGaussTurbulentKinematicViscositySensitivitiesEpsilon *
                       mData.mInvEpsilonSigma;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::VelocityDerivativeData::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
double EpsilonAdjointElementData<TDim, TNumNodes>::ShapeDerivativeData::CalculateEffectiveKinematicViscosityShapeDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::KDerivativeData::CalculateReactionTermDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY

    if (mData.mReactionTerm > 0.0)
    {
        KEpsilonElementData::CalculateGaussGammaTkeSensitivity<TNumNodes>(
            rOutput, mData.mCmu, mData.mGamma, mData.mTurbulentKinematicViscosity,
            mData.mGaussTurbulentKinematicViscositySensitivitiesK, rShapeFunctions);
        noalias(rOutput) = rOutput * mData.mC2;
    }
    else
    {
        rOutput.clear();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::EpsilonDerivativeData::CalculateReactionTermDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY

    if (mData.mReactionTerm > 0.0)
    {
        KEpsilonElementData::CalculateGaussGammaEpsilonSensitivity<TNumNodes>(
            rOutput, mData.mGamma, mData.mTurbulentKinematicViscosity,
            mData.mGaussTurbulentKinematicViscositySensitivitiesEpsilon);
        noalias(rOutput) = rOutput * mData.mC2;
    }
    else
    {
        rOutput.clear();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::VelocityDerivativeData::CalculateReactionTermDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY

    if (mData.mReactionTerm > 0.0)
    {
        noalias(rOutput) = rShapeFunctionDerivatives * (mData.mC1 * 2.0 / 3.0);
    }
    else
    {
        rOutput.clear();
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
double EpsilonAdjointElementData<TDim, TNumNodes>::ShapeDerivativeData::CalculateReactionTermShapeDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    if (mData.mReactionTerm > 0.0)
    {
        return (mData.mC1 * 2.0 / 3.0) * RansCalculationUtilities::GetDivergence(
                                             mData.GetGeometry(), VELOCITY, rDN_Dx_deriv);
    }
    else
    {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::KDerivativeData::CalculateSourceTermDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY

    BoundedVector<double, TNumNodes> production_sensitivities;

    KEpsilonElementData::CalculateProductionScalarSensitivities<TNumNodes>(
        production_sensitivities, mData.mTurbulentKinematicViscosity,
        mData.mProductionTerm, mData.mGaussTurbulentKinematicViscositySensitivitiesK);

    KEpsilonElementData::CalculateGaussGammaTkeSensitivity<TNumNodes>(
        rOutput, mData.mCmu, mData.mGamma, mData.mTurbulentKinematicViscosity,
        mData.mGaussTurbulentKinematicViscositySensitivitiesK, rShapeFunctions);

    noalias(rOutput) =
        rOutput * mData.mProductionTerm + production_sensitivities * mData.mGamma;
    noalias(rOutput) = rOutput * mData.mC1;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::EpsilonDerivativeData::CalculateSourceTermDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::VelocityDerivativeData::CalculateSourceTermDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    KRATOS_TRY

    KEpsilonElementData::CalculateProductionVelocitySensitivities<TDim, TNumNodes>(
        rOutput, mData.mTurbulentKinematicViscosity, mData.mProductionTerm,
        mData.mVelocityGradient, rShapeFunctionDerivatives);

    noalias(rOutput) = rOutput * (mData.mC1 * mData.mGamma);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
double EpsilonAdjointElementData<TDim, TNumNodes>::ShapeDerivativeData::CalculateSourceTermShapeDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    return KEpsilonElementData::CalculateProductionShapeSensitivities<TDim, TNumNodes>(
               mData.mTurbulentKinematicViscosity, 0.0, mData.mProductionTerm,
               mData.mNodalVelocity, rShapeFunctionDerivatives, rDN_Dx_deriv) *
           mData.mC1 * mData.mGamma;
}

// template instantiations
template class EpsilonAdjointElementData<2, 3>;
template class EpsilonAdjointElementData<3, 4>;

} // namespace KEpsilonElementData
///@}
} // namespace Kratos
