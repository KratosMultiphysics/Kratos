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
#include "includes/cfd_variables.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "rans_application_variables.h"

// Derivative data type includes
#include "custom_conditions/data_containers/k_epsilon/epsilon_k_based_wall_condition_data_derivatives.h"

// Include base h
#include "scalar_wall_flux_condition_derivatives.h"

namespace Kratos
{

/***************************************************************************************************/
/***************************************** Static Methods ******************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TConditionDataType>
void ScalarWallFluxConditionDerivatives<TDim, TNumNodes, TConditionDataType>::Check(
    const Condition& rCondition,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    TConditionDataType::Check(rCondition, rProcessInfo);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TConditionDataType>
GeometryData::IntegrationMethod ScalarWallFluxConditionDerivatives<TDim, TNumNodes, TConditionDataType>::GetIntegrationMethod()
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

/***************************************************************************************************/
/********************************************** Data ***********************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TConditionDataType>
ScalarWallFluxConditionDerivatives<TDim, TNumNodes, TConditionDataType>::Data::Data(
    const Condition& rCondition,
    const Element& rParentElement,
    const ProcessInfo& rProcessInfo)
    : mrCondition(rCondition),
      mrParentElement(rParentElement),
      mrConditionData(rCondition.GetGeometry(), rCondition.GetProperties(), rProcessInfo, mParameters)
{
    KRATOS_TRY

    mGaussPointIndex = 0;

    mParameters.mKappa = rProcessInfo[VON_KARMAN];

    const auto& r_element_properties = rParentElement.GetProperties();
    mParameters.mDensity = r_element_properties[DENSITY];
    mParameters.mKinematicViscosity = r_element_properties[DYNAMIC_VISCOSITY] / mParameters.mDensity;

    mrConditionData.CalculateConstants(rProcessInfo);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TConditionDataType>
void ScalarWallFluxConditionDerivatives<TDim, TNumNodes, TConditionDataType>::Data::CalculateGaussPointData(
    const double W,
    const Vector& rN,
    const int Step)
{
    KRATOS_TRY

    mParameters.mYPlus = mrCondition.GetValue(GAUSS_RANS_Y_PLUS)[mGaussPointIndex++];
    mParameters.mWallTurbulentViscosity = mParameters.mKappa * mParameters.mYPlus * mParameters.mKinematicViscosity;

    mrConditionData.CalculateGaussPointData(rN, Step);

    KRATOS_CATCH("");
}

/***************************************************************************************************/
/************************************* Residual Contributions **************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TConditionDataType>
ScalarWallFluxConditionDerivatives<TDim, TNumNodes, TConditionDataType>::ResidualsContributions::ResidualsContributions(
    Data& rData)
    : mrData(rData)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TConditionDataType>
void ScalarWallFluxConditionDerivatives<TDim, TNumNodes, TConditionDataType>::ResidualsContributions::AddGaussPointResidualsContributions(
    VectorN& rResidual,
    const double W,
    const Vector& rN)
{
    noalias(rResidual) += rN * (W * mrData.mrConditionData.GetWallFlux());
}

/***************************************************************************************************/
/************************************** Variable Derivatives ***************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TConditionDataType>
template <class TDerivativesType>
ScalarWallFluxConditionDerivatives<TDim, TNumNodes, TConditionDataType>::VariableDerivatives<
    TDerivativesType>::VariableDerivatives(Data& rData)
    : mrData(rData),
      mResidualWeightDerivativeContributions(rData)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TConditionDataType>
template <class TDerivativesType>
void ScalarWallFluxConditionDerivatives<TDim, TNumNodes, TConditionDataType>::VariableDerivatives<TDerivativesType>::CalculateGaussPointResidualsDerivativeContributions(
    VectorN& rResidualDerivative,
    const IndexType ConditionNodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const double WDerivative,
    const double DetJDerivative,
    const IndexType ParentElementNodeIndex)
{
    KRATOS_TRY

    const TDerivativesType derivative(mrData.mrConditionData);

    const double flux_derivative = derivative.CalculateWallFluxDerivative(
        ConditionNodeIndex, DirectionIndex, W, rN, WDerivative, DetJDerivative,
        ParentElementNodeIndex);

    noalias(rResidualDerivative) = rN * (W * flux_derivative);

    mResidualWeightDerivativeContributions.AddGaussPointResidualsContributions(rResidualDerivative, WDerivative, rN);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TConditionDataType>
template <class TDerivativesType>
void ScalarWallFluxConditionDerivatives<TDim, TNumNodes, TConditionDataType>::VariableDerivatives<TDerivativesType>::CalculateGaussPointResidualsDerivativeContributions(
    VectorN& rResidualDerivative,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const IndexType ParentElementNodeIndex)
{
    KRATOS_TRY

    const TDerivativesType derivative(mrData.mrConditionData);

    const double flux_derivative = derivative.CalculateWallFluxDerivative(
        DirectionIndex, W, rN, ParentElementNodeIndex);

    noalias(rResidualDerivative) = rN * (W * flux_derivative);

    KRATOS_CATCH("");
}

// template instantiations

template class ScalarWallFluxConditionDerivatives<2, 2, typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<2>::Data>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<2>::UDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<2>::KDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<2>::EpsilonDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<2>::ShapeDerivative>;

template class ScalarWallFluxConditionDerivatives<3, 3, typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<3>::Data>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<3>::UDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<3>::KDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<3>::EpsilonDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonKBasedWallConditionDataDerivatives<3>::ShapeDerivative>;

} // namespace Kratos
