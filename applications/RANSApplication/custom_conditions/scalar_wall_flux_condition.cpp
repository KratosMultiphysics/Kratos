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
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "rans_application_variables.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_conditions/data_containers/k_epsilon/k_vis_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_epsilon/epsilon_k_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_epsilon/epsilon_u_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega/omega_k_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega/omega_u_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega/omega_vis_log_law_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega_sst/omega_k_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega_sst/omega_u_based_wall_condition_data.h"

// Include base h
#include "scalar_wall_flux_condition.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
void ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& CurrentProcessInfo) const
{
    if (rResult.size() != TNumNodes) {
        rResult.resize(TNumNodes, false);
    }

    const Variable<double>& r_variable =
        TScalarWallFluxConditionData::GetScalarVariable();

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rResult[i] = Condition::GetGeometry()[i].GetDof(r_variable).EquationId();
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
void ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& CurrentProcessInfo) const
{
    if (rConditionalDofList.size() != TNumNodes) {
        rConditionalDofList.resize(TNumNodes);
    }

    const Variable<double>& r_variable =
        TScalarWallFluxConditionData::GetScalarVariable();

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rConditionalDofList[i] = Condition::GetGeometry()[i].pGetDof(r_variable);
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
void ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>::GetValuesVector(
    Vector& rValues,
    int Step) const
{
    this->GetFirstDerivativesVector(rValues, Step);
}

template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
void ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>::GetFirstDerivativesVector(
    Vector& rValues,
    int Step) const
{
    if (rValues.size() != TNumNodes) {
        rValues.resize(TNumNodes, false);
    }

    const auto& r_geometry = this->GetGeometry();
    const Variable<double>& r_variable =
        TScalarWallFluxConditionData::GetScalarVariable();

    IndexType local_index = 0;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        rValues[local_index++] =
            r_geometry[i_node].FastGetSolutionStepValue(r_variable, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
void ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>::GetSecondDerivativesVector(
    Vector& rValues,
    int Step) const
{
    if (rValues.size() != TNumNodes) {
        rValues.resize(TNumNodes, false);
    }

    const auto& r_geometry = this->GetGeometry();
    const Variable<double>& r_variable =
        TScalarWallFluxConditionData::GetScalarVariable().GetTimeDerivative();

    IndexType local_index = 0;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        rValues[local_index++] =
            r_geometry[i_node].FastGetSolutionStepValue(r_variable, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
void ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Check sizes and initialize matrix
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes) {
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

    // Calculate RHS
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
void ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }

    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

    if (TScalarWallFluxConditionData::IsWallFluxComputed(this->GetGeometry())) {
        const auto& r_geometry = this->GetGeometry();
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        RansCalculationUtilities::CalculateConditionGeometryData(
            r_geometry, this->GetIntegrationMethod(), gauss_weights, shape_functions);
        const IndexType num_gauss_points = gauss_weights.size();

        TScalarWallFluxConditionData r_current_data(
            r_geometry, this->GetProperties(), rCurrentProcessInfo);
        r_current_data.CalculateConstants(rCurrentProcessInfo);

        const auto& gauss_y_plus = this->GetGeometry().GetValue(GAUSS_RANS_Y_PLUS);

        ScalarWallFluxConditionData::Parameters params;

        params.mKappa = rCurrentProcessInfo[VON_KARMAN];

        const auto& r_element_properties = this->GetGeometry().GetValue(NEIGHBOUR_ELEMENTS)[0].GetProperties();
        params.mDensity = r_element_properties[DENSITY];
        params.mKinematicViscosity = r_element_properties[DYNAMIC_VISCOSITY] / params.mDensity;

        for (IndexType g = 0; g < num_gauss_points; ++g) {
            const Vector& gauss_shape_functions = row(shape_functions, g);

            params.mYPlus = gauss_y_plus[g];
            params.mWallTurbulentViscosity = params.mKappa * params.mYPlus * params.mKinematicViscosity;

            const double flux = r_current_data.CalculateWallFlux(gauss_shape_functions, params);

            noalias(rRightHandSideVector) += gauss_shape_functions * (gauss_weights[g] * flux);
        }
    }

    KRATOS_CATCH("");
}
template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
GeometryData::IntegrationMethod  ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
int ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>::Check(
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int check = BaseType::Check(rCurrentProcessInfo);
    if (TScalarWallFluxConditionData::IsWallFluxComputed(this->GetGeometry())) {
        TScalarWallFluxConditionData::Check(*this, rCurrentProcessInfo);

        KRATOS_ERROR_IF_NOT(this->Has(NEIGHBOUR_ELEMENTS))
            << "NEIGHBOUR_ELEMENTS were not found in condition "
            << this->Info() << ".\n";

        KRATOS_ERROR_IF_NOT(this->GetValue(NEIGHBOUR_ELEMENTS).size() == 1)
            << "More than one parent element was found for condition " << this->Info()
            << " [ number of parents = " << this->GetValue(NEIGHBOUR_ELEMENTS).size()
            << " ].\n";

        KRATOS_ERROR_IF_NOT(this->Has(GAUSS_RANS_Y_PLUS))
            << "GAUSS_RANS_Y_PLUS were not found in condition "
            << this->Info() << ".\n";

        const std::size_t number_of_gauss_points = this->GetGeometry().IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_2);
        KRATOS_ERROR_IF_NOT(this->GetValue(GAUSS_RANS_Y_PLUS).size() == number_of_gauss_points)
            << "GAUSS_RANS_Y_PLUS were not initialized properly in condition "
            << this->Info()
            << " [ GAUSS_RANS_Y_PLUS.size() = " << this->GetValue(GAUSS_RANS_Y_PLUS).size()
            << ", required size = " << number_of_gauss_points << " ].\n";

        KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(VON_KARMAN)) << "VON_KARMAN is not found in process info.\n";

        const auto& r_parent_element = this->GetValue(NEIGHBOUR_ELEMENTS)[0];
        const auto& r_parent_element_properties = r_parent_element.GetProperties();

        KRATOS_ERROR_IF_NOT(r_parent_element_properties.Has(DENSITY))
            << "DENSITY is not found in parent element properties. [ "
            "Properties.Id() = "
            << r_parent_element_properties.Id()
            << ", ParentElement.Id() = " << r_parent_element.Id()
            << ", Condition.Id() = " << this->Id() << " ].\n";

        KRATOS_ERROR_IF_NOT(r_parent_element_properties.Has(DYNAMIC_VISCOSITY))
            << "DYNAMIC_VISCOSITY is not found in parent element properties. [ "
            "Properties.Id() = "
            << r_parent_element_properties.Id()
            << ", ParentElement.Id() = " << r_parent_element.Id()
            << ", Condition.Id() = " << this->Id() << " ].\n";
    }

    return check;

    KRATOS_CATCH("");
}

// template instantiations

template class ScalarWallFluxCondition<2, 2, KEpsilonWallConditionData::KVisBasedWallConditionData>;
template class ScalarWallFluxCondition<2, 2, KEpsilonWallConditionData::EpsilonKBasedWallConditionData>;
template class ScalarWallFluxCondition<2, 2, KEpsilonWallConditionData::EpsilonUBasedWallConditionData>;
template class ScalarWallFluxCondition<3, 3, KEpsilonWallConditionData::KVisBasedWallConditionData>;
template class ScalarWallFluxCondition<3, 3, KEpsilonWallConditionData::EpsilonKBasedWallConditionData>;
template class ScalarWallFluxCondition<3, 3, KEpsilonWallConditionData::EpsilonUBasedWallConditionData>;

template class ScalarWallFluxCondition<2, 2, KOmegaWallConditionData::OmegaKBasedWallConditionData>;
template class ScalarWallFluxCondition<2, 2, KOmegaWallConditionData::OmegaUBasedWallConditionData>;
template class ScalarWallFluxCondition<2, 2, KOmegaWallConditionData::OmegaVisLogBasedWallConditionData>;
template class ScalarWallFluxCondition<3, 3, KOmegaWallConditionData::OmegaKBasedWallConditionData>;
template class ScalarWallFluxCondition<3, 3, KOmegaWallConditionData::OmegaUBasedWallConditionData>;
template class ScalarWallFluxCondition<3, 3, KOmegaWallConditionData::OmegaVisLogBasedWallConditionData>;

template class ScalarWallFluxCondition<2, 2, KOmegaSSTWallConditionData::OmegaKBasedWallConditionData>;
template class ScalarWallFluxCondition<2, 2, KOmegaSSTWallConditionData::OmegaUBasedWallConditionData>;
template class ScalarWallFluxCondition<3, 3, KOmegaSSTWallConditionData::OmegaKBasedWallConditionData>;
template class ScalarWallFluxCondition<3, 3, KOmegaSSTWallConditionData::OmegaUBasedWallConditionData>;

} // namespace Kratos.
