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

// External includes

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_constitutive/fluid_constitutive_law.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"

// Include base h
#include "convection_diffusion_reaction_residual_based_flux_corrected_derivatives.h"

namespace Kratos
{

/***************************************************************************************************/
/***************************************** Static Methods ******************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
int ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::Check(
    const GeometryType& rGeometry,
    const ProcessInfo& rProcessInfo)
{
    return 0;
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::GetIntegrationMethod()
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

/***************************************************************************************************/
/********************************************** Data ***********************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::Data::Data(
    const Element& rElement,
    FluidConstitutiveLaw& rFluidConstitutiveLaw)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::Data::Initialize(
    const ProcessInfo& rProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::Data::CalculateGaussPointData(
    const double W,
    const Vector& rN,
    const Matrix& rdNdX)
{
}

/***************************************************************************************************/
/************************************* Residual Contributions **************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::ResidualContributions::ResidualContributions(
    Data& rData)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::ResidualContributions::Initialize(
    Matrix& rOutput,
    const ProcessInfo& rProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::ResidualContributions::AddResidualContribution(
    BoundedVector<double, TElementLocalSize>& rResidual,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX)
{
}

/***************************************************************************************************/
/************************************** Variable Derivatives ***************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes>
template <class TDerivativesType, unsigned int TEquationOffset, unsigned int TDerivativeOffset>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::
    VariableDerivatives<TDerivativesType, TEquationOffset, TDerivativeOffset>::VariableDerivatives(Data& rData)
    : mrData(rData), mResidualWeightDerivativeContributions(rData)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
template <class TDerivativesType, unsigned int TEquationOffset, unsigned int TDerivativeOffset>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::
    VariableDerivatives<TDerivativesType, TEquationOffset, TDerivativeOffset>::Initialize(
        Matrix& rOutput, const ProcessInfo& rProcessInfo)
{
    mBlockSize = rOutput.size2() / TNumNodes;
    mResidualWeightDerivativeContributions.Initialize(rOutput, rProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
template <class TDerivativesType, unsigned int TEquationOffset, unsigned int TDerivativeOffset>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::
    VariableDerivatives<TDerivativesType, TEquationOffset, TDerivativeOffset>::CalculateResidualDerivative(
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
}

template <unsigned int TDim, unsigned int TNumNodes>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::SecondDerivatives::SecondDerivatives(
    const Element& rElement,
    FluidConstitutiveLaw& rFluidConstitutiveLaw)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::SecondDerivatives::Initialize(
    Matrix& rOutput,
    const ProcessInfo& rProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes>::SecondDerivatives::AddResidualDerivativeContributions(
    Matrix& rOutput,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX)
{
}

} // namespace Kratos
