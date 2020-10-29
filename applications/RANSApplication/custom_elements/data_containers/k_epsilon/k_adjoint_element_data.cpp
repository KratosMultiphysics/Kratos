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
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "k_adjoint_element_data.h"

namespace Kratos
{
namespace KEpsilonAdjointElementData
{

//***********************************************************************************//
//************************** Adjoint Element Data ***********************************//

template <unsigned int TDim, unsigned int TNumNodes>
KAdjointElementData<TDim, TNumNodes>::AdjointElementData::AdjointElementData(
    const GeometryType& rGeometry)
    : BaseType(rGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::AdjointElementData::CalculateGaussPointData(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const int Step)
{
    BaseType::CalculateGaussPointData(rShapeFunctions, rShapeFunctionDerivatives, Step);

    mProductionTerm = this->CalculateSourceTerm(rShapeFunctions, rShapeFunctionDerivatives);
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
void KAdjointElementData<TDim, TNumNodes>::AdjointElementData::Check(
    const GeometryType& rGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    // BaseType::Check(rGeometry, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& KAdjointElementData<TDim, TNumNodes>::AdjointElementData::GetAdjointScalarVariable()
{
    return RANS_SCALAR_1_ADJOINT_1;
}

//************************** Adjoint Element Data ***********************************//
//***********************************************************************************//
//************************** Velocity Derivatives ***********************************//

template <unsigned int TDim, unsigned int TNumNodes>
KAdjointElementData<TDim, TNumNodes>::VelocityDerivatives::VelocityDerivatives(
    const AdjointElementData& rElementData)
    : BaseType(rElementData)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<array_1d<double, 3>>& KAdjointElementData<TDim, TNumNodes>::VelocityDerivatives::GetDerivativeVariable()
{
    return VELOCITY;
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::VelocityDerivatives::CalculateEffectiveVelocityDerivatives(
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
void KAdjointElementData<TDim, TNumNodes>::VelocityDerivatives::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::VelocityDerivatives::CalculateReactionTermDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    if (this->GetElementData().mReactionTerm > 0.0) {
        for (IndexType c = 0; c < TNumNodes; ++c) {
            const auto block_row = c * TDim;
            for (IndexType k = 0; k < TDim; ++k) {
                rOutput[block_row + k] = rShapeFunctionDerivatives(c, k) * (2.0 / 3.0);
            }
        }
    } else {
        rOutput.clear();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::VelocityDerivatives::CalculateSourceTermDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& element_data = this->GetElementData();

    KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionVelocitySensitivities(
        rOutput, element_data.mTurbulentKinematicViscosity, element_data.mProductionTerm,
        element_data.mVelocityGradient, rShapeFunctionDerivatives);
}

//************************** Velocity Derivatives ***********************************//
//***********************************************************************************//
//***************************** K Derivatives ***************************************//

template <unsigned int TDim, unsigned int TNumNodes>
KAdjointElementData<TDim, TNumNodes>::KDerivatives::KDerivatives(
    const AdjointElementData& rElementData)
    : BaseType(rElementData)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& KAdjointElementData<TDim, TNumNodes>::KDerivatives::GetDerivativeVariable()
{
    return TURBULENT_KINETIC_ENERGY;
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::KDerivatives::CalculateEffectiveVelocityDerivatives(
    BoundedMatrix<double, TDerivativesSize, 3>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::KDerivatives::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();
    noalias(rOutput) =
        data.mGaussTurbulentKinematicViscositySensitivitiesK * data.mInvTkeSigma;
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::KDerivatives::CalculateReactionTermDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();

    if (data.mReactionTerm > 0.0) {
        KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateGaussGammaTkeSensitivity(
            rOutput, data.mCmu, data.mGamma, data.mTurbulentKinematicViscosity,
            data.mGaussTurbulentKinematicViscositySensitivitiesK, rShapeFunctions);
    } else {
        rOutput.clear();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::KDerivatives::CalculateSourceTermDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();
    KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionScalarSensitivities(
        rOutput, data.mTurbulentKinematicViscosity, data.mProductionTerm,
        data.mGaussTurbulentKinematicViscositySensitivitiesK);
}

//***************************** K Derivatives ***************************************//
//***********************************************************************************//
//************************** Epsilon Derivatives ************************************//

template <unsigned int TDim, unsigned int TNumNodes>
KAdjointElementData<TDim, TNumNodes>::EpsilonDerivatives::EpsilonDerivatives(
    const AdjointElementData& rElementData)
    : BaseType(rElementData)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& KAdjointElementData<TDim, TNumNodes>::EpsilonDerivatives::GetDerivativeVariable()
{
    return TURBULENT_ENERGY_DISSIPATION_RATE;
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::EpsilonDerivatives::CalculateEffectiveVelocityDerivatives(
    BoundedMatrix<double, TDerivativesSize, 3>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    rOutput.clear();
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::EpsilonDerivatives::CalculateEffectiveKinematicViscosityDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();
    noalias(rOutput) =
        data.mGaussTurbulentKinematicViscositySensitivitiesEpsilon * data.mInvTkeSigma;
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::EpsilonDerivatives::CalculateReactionTermDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();

    if (data.mReactionTerm > 0.0) {
        KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateGaussGammaEpsilonSensitivity(
            rOutput, data.mGamma, data.mTurbulentKinematicViscosity,
            data.mGaussTurbulentKinematicViscositySensitivitiesEpsilon);
    } else {
        rOutput.clear();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::EpsilonDerivatives::CalculateSourceTermDerivatives(
    BoundedVector<double, TDerivativesSize>& rOutput,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const auto& data = this->GetElementData();
    KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionScalarSensitivities(
        rOutput, data.mTurbulentKinematicViscosity, data.mProductionTerm,
        data.mGaussTurbulentKinematicViscositySensitivitiesEpsilon);
}

//************************** Epsilon Derivatives ************************************//
//***********************************************************************************//
//************************ Sensitivity Element Data *********************************//

template <unsigned int TDim, unsigned int TNumNodes>
KAdjointElementData<TDim, TNumNodes>::SensitivityElementData::SensitivityElementData(
    const GeometryType& rGeometry)
    : BaseType(rGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void KAdjointElementData<TDim, TNumNodes>::SensitivityElementData::CalculateGaussPointData(
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
void KAdjointElementData<TDim, TNumNodes>::SensitivityElementData::Check(
    const GeometryType& rGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    // BaseType::Check(rGeometry, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<double>& KAdjointElementData<TDim, TNumNodes>::SensitivityElementData::GetAdjointScalarVariable()
{
    return RANS_SCALAR_1_ADJOINT_1;
}

//************************ Sensitivity Element Data *********************************//
//***********************************************************************************//
//*************************** Shape Sensitivity *************************************//

template <unsigned int TDim, unsigned int TNumNodes>
KAdjointElementData<TDim, TNumNodes>::ShapeSensitivity::ShapeSensitivity(
    const SensitivityElementData& rElementData)
    : BaseType(rElementData)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
const Variable<array_1d<double, 3>>& KAdjointElementData<TDim, TNumNodes>::ShapeSensitivity::GetDerivativeVariable()
{
    return SHAPE_SENSITIVITY;
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, 3> KAdjointElementData<TDim, TNumNodes>::ShapeSensitivity::CalculateEffectiveVelocityDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    return ZeroVector(3);
}

template <unsigned int TDim, unsigned int TNumNodes>
double KAdjointElementData<TDim, TNumNodes>::ShapeSensitivity::CalculateEffectiveKinematicViscosityDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
double KAdjointElementData<TDim, TNumNodes>::ShapeSensitivity::CalculateReactionTermDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    const auto& data = this->GetElementData();

    if (data.mReactionTerm > 0.0) {
        return (2.0 / 3.0) * RansCalculationUtilities::GetDivergence(
                                 data.GetGeometry(), VELOCITY, rDN_Dx_deriv);
    } else {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double KAdjointElementData<TDim, TNumNodes>::ShapeSensitivity::CalculateSourceTermDerivatives(
    const ShapeParameter& rShapeParameters,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const double detJ_deriv,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv) const
{
    const auto& data = this->GetElementData();
    return KEpsilonElementData::AdjointUtilities<TDim, TNumNodes>::CalculateProductionShapeSensitivities(
        data.mTurbulentKinematicViscosity, 0.0, data.mProductionTerm,
        data.mNodalVelocity, rShapeFunctionDerivatives, rDN_Dx_deriv);
}

//*************************** Shape Sensitivity *************************************//
//***********************************************************************************//

// template instantiations
template class KAdjointElementData<2, 3>;
template class KAdjointElementData<3, 4>;

} // namespace KEpsilonAdjointElementData
} // namespace Kratos