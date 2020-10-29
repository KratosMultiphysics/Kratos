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

// Project includes
#include "geometries/geometry_data.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

// Application includes
#include "custom_elements/data_containers/k_epsilon/adjoint_element_data_utilities.h"
#include "custom_elements/data_containers/k_epsilon/element_data_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "epsilon_adjoint_element_data.h"

namespace Kratos
{
namespace KEpsilonAdjointElementData
{

//***********************************************************************************//
//************************** Adjoint Element Data ***********************************//

template <unsigned int TDim, unsigned int TNumNodes>
EpsilonAdjointElementData<TDim, TNumNodes>::AdjointElementData::AdjointElementData(
    const GeometryType& rGeometry)
    : BaseType(rGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::AdjointElementData::CalculateGaussPointData(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const int Step)
{
    BaseType::CalculateGaussPointData(rShapeFunctions, rShapeFunctionDerivatives, Step);

    mProductionTerm = KEpsilonElementData::CalculateSourceTerm<TDim>(this->mVelocityGradient, this->mTurbulentKinematicViscosity);
    mEffectiveKinematicViscosity = this->CalculateEffectiveKinematicViscosity(rShapeFunctions, rShapeFunctionDerivatives);
    mReactionTerm = this->CalculateReactionTerm(rShapeFunctions, rShapeFunctionDerivatives);

    BoundedVector<double, TNumNodes> nodal_tke;
    BoundedVector<double, TNumNodes> nodal_epsilon;

    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
        const auto& r_node = this->GetGeometry()[i_node];
        const array_1d<double, 3>& rVelocity =
            r_node.FastGetSolutionStepValue(VELOCITY, Step);
        for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim) {
            mNodalVelocity(i_node, i_dim) = rVelocity[i_dim];
        }
        nodal_tke[i_node] = r_node.FastGetSolutionStepValue(TURBULENT_KINETIC_ENERGY, Step);
        nodal_epsilon[i_node] = r_node.FastGetSolutionStepValue(TURBULENT_ENERGY_DISSIPATION_RATE, Step);
    }

    KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateNodalNutTKESensitivity(
        this->mGaussTurbulentKinematicViscositySensitivitiesK, this->mCmu,
        nodal_tke, nodal_epsilon);
    KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateNodalNutEpsilonSensitivity(
        this->mGaussTurbulentKinematicViscositySensitivitiesEpsilon, this->mCmu,
        nodal_tke, nodal_epsilon);

    RansCalculationUtilities::CalculateGaussSensitivities<TNumNodes>(
        this->mGaussTurbulentKinematicViscositySensitivitiesK,
        this->mGaussTurbulentKinematicViscositySensitivitiesK, rShapeFunctions);
    RansCalculationUtilities::CalculateGaussSensitivities<TNumNodes>(
        this->mGaussTurbulentKinematicViscositySensitivitiesEpsilon,
        this->mGaussTurbulentKinematicViscositySensitivitiesEpsilon, rShapeFunctions);
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::AdjointElementData::Check(
    const GeometryType& rGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    // BaseType::Check(rGeometry, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& EpsilonAdjointElementData<TDim, TNumNodes>::AdjointElementData::GetAdjointScalarVariable()
{
    return RANS_SCALAR_2_ADJOINT_1;
}

//************************** Adjoint Element Data ***********************************//
//***********************************************************************************//
//************************** Velocity Derivatives ***********************************//

template <unsigned int TDim, unsigned int TNumNodes>
EpsilonAdjointElementData<TDim, TNumNodes>::VelocityDerivatives::VelocityDerivatives(
    const AdjointElementData& rElementData)
    : BaseType(rElementData)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<array_1d<double, 3>>& EpsilonAdjointElementData<TDim, TNumNodes>::VelocityDerivatives::GetDerivativeVariable()
{
    return VELOCITY;
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::VelocityDerivatives::CalculateEffectiveVelocityDerivatives(
    BoundedMatrix<double, TDerivativesSize, 3>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();

    for (IndexType c = 0; c < TNumNodes; ++c) {
        const auto block_row = c * TDim;
        for (IndexType k = 0; k < TDim; ++k) {
            rOutput(block_row + k, k) = rShapeFunctions[c];
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::VelocityDerivatives::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::VelocityDerivatives::CalculateReactionTermDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();
    if (data.mReactionTerm > 0.0) {
        for (IndexType c = 0; c < TNumNodes; ++c) {
            const auto block_row = c * TDim;
            for (IndexType k = 0; k < TDim; ++k) {
                rOutput[block_row + k] = data.mC1 * rShapeFunctionDerivatives(c, k) * (2.0 / 3.0);
            }
        }
    } else {
        rOutput.clear();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::VelocityDerivatives::CalculateSourceTermDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();

    KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionVelocitySensitivities(
        rOutput, data.mTurbulentKinematicViscosity, data.mProductionTerm,
        data.mVelocityGradient, rShapeFunctionDerivatives);
    noalias(rOutput) = rOutput * (data.mC1 * data.mGamma);
}

//************************** Velocity Derivatives ***********************************//
//***********************************************************************************//
//***************************** K Derivatives ***************************************//

template <unsigned int TDim, unsigned int TNumNodes>
EpsilonAdjointElementData<TDim, TNumNodes>::KDerivatives::KDerivatives(
    const AdjointElementData& rElementData)
    : BaseType(rElementData)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& EpsilonAdjointElementData<TDim, TNumNodes>::KDerivatives::GetDerivativeVariable()
{
    return TURBULENT_KINETIC_ENERGY;
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::KDerivatives::CalculateEffectiveVelocityDerivatives(
    BoundedMatrix<double, TDerivativesSize, 3>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::KDerivatives::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();
    noalias(rOutput) =
        data.mGaussTurbulentKinematicViscositySensitivitiesK * data.mInvEpsilonSigma;
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::KDerivatives::CalculateReactionTermDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();

    if (data.mReactionTerm > 0.0) {
        KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateGaussGammaTkeSensitivity(
            rOutput, data.mCmu, data.mGamma, data.mTurbulentKinematicViscosity,
            data.mGaussTurbulentKinematicViscositySensitivitiesK, rShapeFunctions);
        noalias(rOutput) = rOutput * data.mC2;
    } else {
        rOutput.clear();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::KDerivatives::CalculateSourceTermDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();

    BoundedVector<double, TDerivativesSize> production_sensitivities;

    KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionScalarSensitivities(
        production_sensitivities, data.mTurbulentKinematicViscosity, data.mProductionTerm,
        data.mGaussTurbulentKinematicViscositySensitivitiesK);

    KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateGaussGammaTkeSensitivity(
        rOutput, data.mCmu, data.mGamma, data.mTurbulentKinematicViscosity,
        data.mGaussTurbulentKinematicViscositySensitivitiesK, rShapeFunctions);

    noalias(rOutput) = rOutput * data.mProductionTerm + production_sensitivities * data.mGamma;
    noalias(rOutput) = rOutput * data.mC1;
}

//***************************** K Derivatives ***************************************//
//***********************************************************************************//
//************************** Epsilon Derivatives ************************************//

template <unsigned int TDim, unsigned int TNumNodes>
EpsilonAdjointElementData<TDim, TNumNodes>::EpsilonDerivatives::EpsilonDerivatives(
    const AdjointElementData& rElementData)
    : BaseType(rElementData)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& EpsilonAdjointElementData<TDim, TNumNodes>::EpsilonDerivatives::GetDerivativeVariable()
{
    return TURBULENT_ENERGY_DISSIPATION_RATE;
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::EpsilonDerivatives::CalculateEffectiveVelocityDerivatives(
    BoundedMatrix<double, TDerivativesSize, 3>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::EpsilonDerivatives::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();
    noalias(rOutput) =
        data.mGaussTurbulentKinematicViscositySensitivitiesEpsilon * data.mInvEpsilonSigma;
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::EpsilonDerivatives::CalculateReactionTermDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();

    if (data.mReactionTerm > 0.0) {
        KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateGaussGammaEpsilonSensitivity(
            rOutput, data.mGamma, data.mTurbulentKinematicViscosity,
            data.mGaussTurbulentKinematicViscositySensitivitiesEpsilon);
        noalias(rOutput) = rOutput * data.mC2;
    } else {
        rOutput.clear();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::EpsilonDerivatives::CalculateSourceTermDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();
}

//************************** Epsilon Derivatives ************************************//
//***********************************************************************************//
//************************ Sensitivity Element Data *********************************//

template <unsigned int TDim, unsigned int TNumNodes>
EpsilonAdjointElementData<TDim, TNumNodes>::SensitivityElementData::SensitivityElementData(
    const GeometryType& rGeometry)
    : BaseType(rGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::SensitivityElementData::CalculateGaussPointData(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const int Step)
{
    BaseType::CalculateGaussPointData(rShapeFunctions, rShapeFunctionDerivatives, Step);

    mProductionTerm = this->CalculateSourceTerm(rShapeFunctions, rShapeFunctionDerivatives);
    mReactionTerm = this->CalculateReactionTerm(rShapeFunctions, rShapeFunctionDerivatives);

    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
        const auto& r_node = this->GetGeometry()[i_node];
        const array_1d<double, 3>& rVelocity =
            r_node.FastGetSolutionStepValue(VELOCITY, Step);
        for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim) {
            mNodalVelocity(i_node, i_dim) = rVelocity[i_dim];
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void EpsilonAdjointElementData<TDim, TNumNodes>::SensitivityElementData::Check(
    const GeometryType& rGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    // BaseType::Check(rGeometry, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& EpsilonAdjointElementData<TDim, TNumNodes>::SensitivityElementData::GetAdjointScalarVariable()
{
    return RANS_SCALAR_2_ADJOINT_1;
}

//************************ Sensitivity Element Data *********************************//
//***********************************************************************************//
//*************************** Shape Sensitivity *************************************//

template <unsigned int TDim, unsigned int TNumNodes>
EpsilonAdjointElementData<TDim, TNumNodes>::ShapeSensitivity::ShapeSensitivity(
    const SensitivityElementData& rElementData)
    : BaseType(rElementData)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<array_1d<double, 3>>& EpsilonAdjointElementData<TDim, TNumNodes>::ShapeSensitivity::GetDerivativeVariable()
{
    return SHAPE_SENSITIVITY;
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, 3> EpsilonAdjointElementData<TDim, TNumNodes>::ShapeSensitivity::CalculateEffectiveVelocityDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    return ZeroVector(3);
}

template <unsigned int TDim, unsigned int TNumNodes>
double EpsilonAdjointElementData<TDim, TNumNodes>::ShapeSensitivity::CalculateEffectiveKinematicViscosityDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
double EpsilonAdjointElementData<TDim, TNumNodes>::ShapeSensitivity::CalculateReactionTermDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    const auto& data = this->GetElementData();

    if (data.mReactionTerm > 0.0) {
        return (data.mC1 * 2.0 / 3.0) * RansCalculationUtilities::GetDivergence(
                                             data.GetGeometry(), VELOCITY, rDN_Dx_deriv);
    } else {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double EpsilonAdjointElementData<TDim, TNumNodes>::ShapeSensitivity::CalculateSourceTermDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    const auto& data = this->GetElementData();
    return KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionShapeSensitivities(
               data.mTurbulentKinematicViscosity, 0.0, data.mProductionTerm,
               data.mNodalVelocity, rShapeFunctionDerivatives, rDN_Dx_deriv) *
           data.mC1 * data.mGamma;
}

//*************************** Shape Sensitivity *************************************//
//***********************************************************************************//



// template instantiations
template class EpsilonAdjointElementData<2, 3>;
template class EpsilonAdjointElementData<3, 4>;

} // namespace KEpsilonAdjointElementData
} // namespace Kratos