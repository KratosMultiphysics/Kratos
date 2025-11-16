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
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "rans_application_variables.h"

// Derivative data type includes
// k-epsilon
#include "custom_conditions/data_containers/k_epsilon/epsilon_k_based_wall_condition_data_derivatives.h"
#include "custom_conditions/data_containers/k_epsilon/epsilon_u_based_wall_condition_data_derivatives.h"

// k-omega
#include "custom_conditions/data_containers/k_omega/omega_k_based_wall_condition_data_derivatives.h"
#include "custom_conditions/data_containers/k_omega/omega_u_based_wall_condition_data_derivatives.h"

// k-omega-sst
#include "custom_conditions/data_containers/k_omega_sst/omega_k_based_wall_condition_data_derivatives.h"
#include "custom_conditions/data_containers/k_omega_sst/omega_u_based_wall_condition_data_derivatives.h"

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

    KRATOS_ERROR_IF_NOT(rCondition.Has(NEIGHBOUR_ELEMENTS))
        << "NEIGHBOUR_ELEMENTS were not found in condition "
        << rCondition.Info() << ".\n";

    KRATOS_ERROR_IF_NOT(rCondition.GetValue(NEIGHBOUR_ELEMENTS).size() == 1)
        << "More than one parent element was found for condition " << rCondition.Info()
        << " [ number of parents = " << rCondition.GetValue(NEIGHBOUR_ELEMENTS).size()
        << " ].\n";

    KRATOS_ERROR_IF_NOT(rCondition.Has(GAUSS_RANS_Y_PLUS))
        << "GAUSS_RANS_Y_PLUS were not found in condition "
        << rCondition.Info() << ".\n";

    const IndexType number_of_gauss_points = rCondition.GetGeometry().IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_2);
    KRATOS_ERROR_IF_NOT(rCondition.GetValue(GAUSS_RANS_Y_PLUS).size() == number_of_gauss_points)
        << "GAUSS_RANS_Y_PLUS were not initialized properly in condition "
        << rCondition.Info()
        << " [ GAUSS_RANS_Y_PLUS.size() = " << rCondition.GetValue(GAUSS_RANS_Y_PLUS).size()
        << ", required size = " << number_of_gauss_points << " ].\n";

    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(VON_KARMAN)) << "VON_KARMAN is not found in process info.\n";

    const auto& r_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];
    const auto& r_parent_element_properties = r_parent_element.GetProperties();

    KRATOS_ERROR_IF_NOT(r_parent_element_properties.Has(DENSITY))
        << "DENSITY is not found in parent element properties. [ "
        "Properties.Id() = "
        << r_parent_element_properties.Id()
        << ", ParentElement.Id() = " << r_parent_element.Id()
        << ", Condition.Id() = " << rCondition.Id() << " ].\n";

    KRATOS_ERROR_IF_NOT(r_parent_element_properties.Has(DYNAMIC_VISCOSITY))
        << "DYNAMIC_VISCOSITY is not found in parent element properties. [ "
        "Properties.Id() = "
        << r_parent_element_properties.Id()
        << ", ParentElement.Id() = " << r_parent_element.Id()
        << ", Condition.Id() = " << rCondition.Id() << " ].\n";

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TConditionDataType>
GeometryData::IntegrationMethod ScalarWallFluxConditionDerivatives<TDim, TNumNodes, TConditionDataType>::GetIntegrationMethod()
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}
template <unsigned int TDim, unsigned int TNumNodes, class TConditionDataType>
void ScalarWallFluxConditionDerivatives<TDim, TNumNodes, TConditionDataType>::InitializeCondition(
    Condition& rCondition,
    const ProcessInfo& rProcessInfo)
{
    TConditionDataType::InitializeCondition(rCondition, rProcessInfo);
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
    const IndexType ParentElementNodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN,
    const double WDerivative,
    const double DetJDerivative)
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
    const IndexType ParentElementNodeIndex,
    const IndexType DirectionIndex,
    const double W,
    const Vector& rN)
{
    KRATOS_TRY

    const TDerivativesType derivative(mrData.mrConditionData);

    const double flux_derivative = derivative.CalculateWallFluxDerivative(
        DirectionIndex, W, rN, ParentElementNodeIndex);

    noalias(rResidualDerivative) = rN * (W * flux_derivative);

    KRATOS_CATCH("");
}

// template instantiations

// k-epsilon k based
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

// k-epsilon u based
template class ScalarWallFluxConditionDerivatives<2, 2, typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<2>::Data>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<2>::UDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<2>::KDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<2>::EpsilonDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<2>::ShapeDerivative>;

template class ScalarWallFluxConditionDerivatives<3, 3, typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<3>::Data>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<3>::UDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<3>::KDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<3>::EpsilonDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KEpsilonWallConditionData::EpsilonUBasedWallConditionDataDerivatives<3>::ShapeDerivative>;

// k-omega k based
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::Data>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::UDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::KDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::OmegaDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::ShapeDerivative>;

template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::Data>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::UDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::KDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::OmegaDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::ShapeDerivative>;

// k-omega u based
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::Data>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::UDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::KDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::OmegaDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::ShapeDerivative>;

template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::Data>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::UDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::KDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::OmegaDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::ShapeDerivative>;

// k-omega-sst k based
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::Data>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::UDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::KDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::OmegaDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<2>::ShapeDerivative>;

template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::Data>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::UDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::KDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::OmegaDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaKBasedWallConditionDataDerivatives<3>::ShapeDerivative>;

// k-omega-sst u based
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::Data>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::UDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::KDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::OmegaDerivative>;
template class ScalarWallFluxConditionDerivatives<2, 2, typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<2>::ShapeDerivative>;

template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::Data>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::UDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::KDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::OmegaDerivative>;
template class ScalarWallFluxConditionDerivatives<3, 3, typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::Data>::VariableDerivatives<typename KOmegaSSTWallConditionData::OmegaUBasedWallConditionDataDerivatives<3>::ShapeDerivative>;

} // namespace Kratos
