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

// Derivative data type includes
#include "custom_elements/data_containers/k_epsilon/k_element_data_derivatives.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data_derivatives.h"

// Include base h
#include "convection_diffusion_reaction_residual_based_flux_corrected_derivatives.h"

namespace Kratos
{

/***************************************************************************************************/
/***************************************** Static Methods ******************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
int ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::Check(
    const GeometryType& rGeometry,
    const ProcessInfo& rProcessInfo)
{
    return 0;
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
GeometryData::IntegrationMethod ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::GetIntegrationMethod()
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

/***************************************************************************************************/
/********************************************** Data ***********************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::Data::Data(
    const Element& rElement,
    const ProcessInfo& rProcessInfo)
    : mrElement(rElement),
      mrElementData(rElement.GetGeometry(), rElement.GetProperties(), rProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::Data::Initialize(
    const ProcessInfo& rProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::Data::CalculateGaussPointData(
    const double W,
    const Vector& rN,
    const Matrix& rdNdX)
{
}

/***************************************************************************************************/
/************************************* Residual Contributions **************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::ResidualContributions::ResidualContributions(
    Data& rData)
    : mrData(rData)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::ResidualContributions::Initialize(
    Vector& rResidual,
    const ProcessInfo& rProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::ResidualContributions::AddResidualContribution(
    Vector& rResidual,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::ResidualContributions::Finalize(
    Vector& rResidual,
    const ProcessInfo& rProcessInfo)
{
}

/***************************************************************************************************/
/**************************************** First Derivative *****************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
template <class TDerivativesType>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::FirstDerivative<
    TDerivativesType>::FirstDerivative(Data& rData)
    : mrData(rData),
      mResidualWeightDerivativeContributions(rData)
{
    // clear the matrices holding primal damping matrix derivatives
    for (IndexType i = 0; i < TDerivativesSize; ++i) {
        mPrimalMatrixDerivatives[i].clear();
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
template <class TDerivativesType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::FirstDerivative<TDerivativesType>::Initialize(
    Vector& rResidualDerivative,
    const ProcessInfo& rProcessInfo)
{
    mBlockSize = rResidualDerivative.size() / TNumNodes;
    mResidualWeightDerivativeContributions.Initialize(rResidualDerivative, rProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
template <class TDerivativesType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::FirstDerivative<TDerivativesType>::CalculateResidualDerivative(
    Vector& rResidualDerivative,
    const int NodeIndex,
    const int DirectionIndex,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX,
    const double WDerivative,
    const double DetJDerivative,
    const Matrix& rdNdXDerivative)
{
    KRATOS_TRY

    const auto& element_data = mrData.mrElementData;

    // create derivatives object
    const TDerivativesType derivative(element_data, NodeIndex, DirectionIndex, W, rN,
                                      rdNdX, WDerivative, DetJDerivative, rdNdXDerivative);

    // compute primal equation coefficients
    const auto& velocity = element_data.GetEffectiveVelocity();
    const auto& viscosity = element_data.GetEffectiveKinematicViscosity();
    const auto& reaction_term = element_data.GetReactionTerm();
    const auto& source_term = element_data.GetSourceTerm();

    // compute primal equation coefficient derivatives
    const auto& velocity_derivative = derivative.CalculateEffectiveVelocityDerivative();
    const auto& viscosity_derivative = derivative.CalculateEffectiveKinematicViscosityDerivative();
    const auto& reaction_term_derivative = derivative.CalculateReactionTermDerivative();
    const auto& source_term_derivative = derivative.CalculateSourceTermDerivative();

    MatrixNN& r_primal_matrix_derivative = mPrimalMatrixDerivatives[TDerivativeDimension * NodeIndex + DirectionIndex];

    for (IndexType a = 0; a < TNumNodes; ++a) {
        // compute LHS matrix derivative for residual based flux corrected stabilization
        for (IndexType b = 0; b < TNumNodes; ++b) {

        }

        // compute the residual derivative without discrete upwind operator and positivity preserving
        // coefficient contributions to derivatives.
        double value = 0.0;



        rResidualDerivative[a] = value;
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
template <class TDerivativesType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::FirstDerivative<TDerivativesType>::Finalize(
    Vector& rResidualDerivative,
    const ProcessInfo& rProcessInfo)
{
    mResidualWeightDerivativeContributions.Finalize(rResidualDerivative, rProcessInfo);
}

/***************************************************************************************************/
/*************************************** Second Derivative *****************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::SecondDerivative::SecondDerivative(
    Data& rData)
    : mrData(rData)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::SecondDerivative::Initialize(
    Vector& rResidualDerivative,
    const ProcessInfo& rProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::SecondDerivative::AddResidualDerivativeContributions(
    Vector& rResidualDerivative,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::SecondDerivative::Finalize(
    Vector& rResidualDerivative,
    const ProcessInfo& rProcessInfo)
{
}

// template instantiations

// k-epsilon k element derivatives
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::KElementDataDerivatives<2, 3>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::KElementDataDerivatives<2, 3>::Data>::FirstDerivative<KEpsilonElementData::KElementDataDerivatives<2, 3>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::KElementDataDerivatives<2, 3>::Data>::FirstDerivative<KEpsilonElementData::KElementDataDerivatives<2, 3>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::KElementDataDerivatives<2, 3>::Data>::FirstDerivative<KEpsilonElementData::KElementDataDerivatives<2, 3>::EpsilonDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::KElementDataDerivatives<2, 3>::Data>::FirstDerivative<KEpsilonElementData::KElementDataDerivatives<2, 3>::ShapeDerivative>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::KElementDataDerivatives<3, 4>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::KElementDataDerivatives<3, 4>::Data>::FirstDerivative<KEpsilonElementData::KElementDataDerivatives<3, 4>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::KElementDataDerivatives<3, 4>::Data>::FirstDerivative<KEpsilonElementData::KElementDataDerivatives<3, 4>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::KElementDataDerivatives<3, 4>::Data>::FirstDerivative<KEpsilonElementData::KElementDataDerivatives<3, 4>::EpsilonDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::KElementDataDerivatives<3, 4>::Data>::FirstDerivative<KEpsilonElementData::KElementDataDerivatives<3, 4>::ShapeDerivative>;

// k-epsilon epsilon element derivatives
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::Data>::FirstDerivative<KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::Data>::FirstDerivative<KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::Data>::FirstDerivative<KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::EpsilonDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::Data>::FirstDerivative<KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::ShapeDerivative>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::Data>::FirstDerivative<KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::Data>::FirstDerivative<KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::Data>::FirstDerivative<KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::EpsilonDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::Data>::FirstDerivative<KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::ShapeDerivative>;

} // namespace Kratos
