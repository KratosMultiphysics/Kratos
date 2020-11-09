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
#include "includes/checks.h"
#include "includes/define.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_conditions/data_containers/k_epsilon/epsilon_k_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_epsilon/epsilon_u_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega/omega_k_based_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega/omega_u_based_wall_condition_data.h"

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

    if (RansCalculationUtilities::IsWallFunctionActive(*this)) {
        const auto& r_geometry = this->GetGeometry();
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        RansCalculationUtilities::CalculateConditionGeometryData(
            r_geometry, TScalarWallFluxConditionData::GetIntegrationMethod(),
            gauss_weights, shape_functions);
        const IndexType num_gauss_points = gauss_weights.size();

        TScalarWallFluxConditionData r_current_data(r_geometry);

        r_current_data.CalculateConstants(rCurrentProcessInfo);

        if (r_current_data.IsWallFluxComputable()) {
            for (IndexType g = 0; g < num_gauss_points; ++g) {
                const Vector& gauss_shape_functions = row(shape_functions, g);

                const double flux = r_current_data.CalculateWallFlux(gauss_shape_functions);

                noalias(rRightHandSideVector) +=
                    gauss_shape_functions * (gauss_weights[g] * flux);
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TScalarWallFluxConditionData>
int ScalarWallFluxCondition<TDim, TNumNodes, TScalarWallFluxConditionData>::Check(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int check = BaseType::Check(rCurrentProcessInfo);
    TScalarWallFluxConditionData::Check(this->GetGeometry(), rCurrentProcessInfo);

    return check;

    KRATOS_CATCH("");
}

// template instantiations

template class ScalarWallFluxCondition<2, 2, KEpsilonWallConditionData::EpsilonKBasedWallConditionData>;
template class ScalarWallFluxCondition<2, 2, KEpsilonWallConditionData::EpsilonUBasedWallConditionData>;
template class ScalarWallFluxCondition<3, 3, KEpsilonWallConditionData::EpsilonKBasedWallConditionData>;
template class ScalarWallFluxCondition<3, 3, KEpsilonWallConditionData::EpsilonUBasedWallConditionData>;

template class ScalarWallFluxCondition<2, 2, KOmegaWallConditionData::OmegaKBasedWallConditionData>;
template class ScalarWallFluxCondition<2, 2, KOmegaWallConditionData::OmegaUBasedWallConditionData>;
template class ScalarWallFluxCondition<3, 3, KOmegaWallConditionData::OmegaKBasedWallConditionData>;
template class ScalarWallFluxCondition<3, 3, KOmegaWallConditionData::OmegaUBasedWallConditionData>;

} // namespace Kratos.
