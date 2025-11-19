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
#include "includes/element.h"
#include "includes/properties.h"
#include "utilities/adjoint_extensions.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/indirect_scalar.h"

// Application includes
#include "custom_utilities/rans_adjoint_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "fluid_dynamics_application_variables.h"
#include "rans_application_variables.h"

// Data container includes
// k-epsilon
#include "custom_conditions/data_containers/k_epsilon/vms_monolithic_kk_based_epsilon_wall_condition_data.h"
#include "custom_conditions/data_containers/k_epsilon/vms_monolithic_ku_based_epsilon_wall_condition_data.h"

// k-omega
#include "custom_conditions/data_containers/k_omega/vms_monolithic_kk_based_omega_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega/vms_monolithic_ku_based_omega_wall_condition_data.h"

// k-omega-sst
#include "custom_conditions/data_containers/k_omega_sst/vms_monolithic_kk_based_omega_wall_condition_data.h"
#include "custom_conditions/data_containers/k_omega_sst/vms_monolithic_ku_based_omega_wall_condition_data.h"

// Include base h
#include "two_equation_turbulence_model_adjoint_condition.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::ThisExtensions::ThisExtensions(Condition* pCondition)
    : mpCondition{pCondition}
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::ThisExtensions::GetFirstDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpCondition->GetGeometry()[NodeId];
    rVector.resize(TDim + 3);

    const auto& dofs_list = TAdjointConditionData::GetDofVariablesList();

    for (unsigned int i = 0; i < TDim; ++i) {
        rVector[i] = MakeIndirectScalar(r_node, (*dofs_list[i]).GetTimeDerivative(), Step);
    }

    rVector[TDim] = IndirectScalar<double>{}; // pressure
    rVector[TDim + 1] = MakeIndirectScalar(r_node, (*dofs_list[TDim + 1]).GetTimeDerivative(), Step);
    rVector[TDim + 2] = MakeIndirectScalar(r_node, (*dofs_list[TDim + 2]).GetTimeDerivative(), Step);

}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::ThisExtensions::GetSecondDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpCondition->GetGeometry()[NodeId];
    rVector.resize(TDim + 3);

    const auto& dofs_list = TAdjointConditionData::GetDofVariablesList();

    for (unsigned int i = 0; i < TDim; ++i) {
        rVector[i] = MakeIndirectScalar(r_node, (*dofs_list[i]).GetTimeDerivative().GetTimeDerivative(), Step);
    }

    rVector[TDim] = IndirectScalar<double>{}; // pressure
    rVector[TDim + 1] = MakeIndirectScalar(r_node, (*dofs_list[TDim + 1]).GetTimeDerivative().GetTimeDerivative(), Step);
    rVector[TDim + 2] = MakeIndirectScalar(r_node, (*dofs_list[TDim + 2]).GetTimeDerivative().GetTimeDerivative(), Step);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::ThisExtensions::GetAuxiliaryVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpCondition->GetGeometry()[NodeId];
    rVector.resize(TDim + 3);

    const auto& dofs_list = TAdjointConditionData::GetDofVariablesList();

    for (unsigned int i = 0; i < TDim; ++i) {
        rVector[i] = MakeIndirectScalar(r_node, (*dofs_list[i]).GetTimeDerivative().GetTimeDerivative().GetTimeDerivative(), Step);
    }

    rVector[TDim] = IndirectScalar<double>{}; // pressure
    rVector[TDim + 1] = MakeIndirectScalar(r_node, (*dofs_list[TDim + 1]).GetTimeDerivative().GetTimeDerivative().GetTimeDerivative(), Step);
    rVector[TDim + 2] = MakeIndirectScalar(r_node, (*dofs_list[TDim + 2]).GetTimeDerivative().GetTimeDerivative().GetTimeDerivative(), Step);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::ThisExtensions::GetFirstDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(3);
    rVariables[0] = &ADJOINT_FLUID_VECTOR_2;
    rVariables[1] = &RANS_SCALAR_1_ADJOINT_2;
    rVariables[2] = &RANS_SCALAR_2_ADJOINT_2;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::ThisExtensions::GetSecondDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(3);
    rVariables[0] = &ADJOINT_FLUID_VECTOR_3;
    rVariables[1] = &RANS_SCALAR_1_ADJOINT_3;
    rVariables[2] = &RANS_SCALAR_2_ADJOINT_3;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::ThisExtensions::GetAuxiliaryVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(3);
    rVariables[0] = &AUX_ADJOINT_FLUID_VECTOR_1;
    rVariables[1] = &RANS_AUX_ADJOINT_SCALAR_1;
    rVariables[2] = &RANS_AUX_ADJOINT_SCALAR_2;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::TwoEquationTurbulenceModelAdjointCondition(IndexType NewId)
    : BaseType(NewId)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::TwoEquationTurbulenceModelAdjointCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::TwoEquationTurbulenceModelAdjointCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::~TwoEquationTurbulenceModelAdjointCondition()
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
Condition::Pointer TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>>(
        NewId, Condition::GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
Condition::Pointer TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
Condition::Pointer TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::Clone(
    IndexType NewId,
    NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>>(
        NewId, Condition::GetGeometry().Create(ThisNodes), Condition::pGetProperties());
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
int TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
     if (RansCalculationUtilities::IsWallFunctionActive(this->GetGeometry())) {
        TAdjointConditionData::Check(*this, rCurrentProcessInfo);
     }
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::EquationIdVector(
    EquationIdVectorType& rConditionalEquationIdList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rConditionalEquationIdList.size() != TConditionLocalSize) {
        rConditionalEquationIdList.resize(TConditionLocalSize, false);
    }

    const auto& r_variables_list = TAdjointConditionData::GetDofVariablesList();

    IndexType local_index = 0;
    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        for (const auto p_variable : r_variables_list) {
            rConditionalEquationIdList[local_index++] = r_node.GetDof(*p_variable).EquationId();
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_variables_list = TAdjointConditionData::GetDofVariablesList();

    if (rCurrentProcessInfo[FRACTIONAL_STEP] == 200 && RansCalculationUtilities::IsWallFunctionActive(this->GetGeometry())) {
        // add dofs for the shape sensitivity calculation
        // it is required to add all the nodes corresponding to parent element
        // because, wall distance calculation involves parent nodes.
        const auto& r_parent_element_geometry = this->GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry();
        const IndexType coords_size = (TDim + 3) * r_parent_element_geometry.size();

        BoundedVector<IndexType, TNumNodes> condition_node_parent_index;
        std::vector<IndexType> parent_only_node_indices;

        ComputeParentElementNodesToConditionNodesMap(
            condition_node_parent_index, parent_only_node_indices,
            this->GetGeometry(), r_parent_element_geometry);


        if (rConditionalDofList.size() != coords_size) {
            rConditionalDofList.resize(coords_size);
        }

        IndexType local_index = 0;

        // add condition dofs
        for (IndexType i = 0; i < TNumNodes; ++i) {
            const auto& r_node = this->GetGeometry()[i];
            for (const auto p_variable : r_variables_list) {
                rConditionalDofList[local_index++] = r_node.pGetDof(*p_variable);
            }
        }

        // add parent element dofs
        for (IndexType i = 0; i < parent_only_node_indices.size(); ++i) {
            const auto& r_node = r_parent_element_geometry[parent_only_node_indices[i]];
            for (const auto p_variable : r_variables_list) {
                rConditionalDofList[local_index++] = r_node.pGetDof(*p_variable);
            }
        }
    } else {
        // add dofs for adjoint solve
        if (rConditionalDofList.size() != TConditionLocalSize) {
            rConditionalDofList.resize(TConditionLocalSize);
        }

        IndexType local_index = 0;
        for (IndexType i = 0; i < TNumNodes; ++i) {
            const auto& r_node = this->GetGeometry()[i];
            for (const auto p_variable : r_variables_list) {
                rConditionalDofList[local_index++] = r_node.pGetDof(*p_variable);
            }
        }
    }


}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::GetValuesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TConditionLocalSize) {
        rValues.resize(TConditionLocalSize, false);
    }

    const auto& r_variables_list = TAdjointConditionData::GetDofVariablesList();

    IndexType local_index = 0;
    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        for (const auto p_variable : r_variables_list) {
            rValues[local_index++] = r_node.FastGetSolutionStepValue(*p_variable, Step);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::GetFirstDerivativesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TConditionLocalSize) {
        rValues.resize(TConditionLocalSize, false);
    }
    rValues.clear();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::GetSecondDerivativesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TConditionLocalSize) {
        rValues.resize(TConditionLocalSize, false);
    }

    const auto& r_variables_list = TAdjointConditionData::GetDofVariablesList();

    IndexType local_index = 0;
    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        for (IndexType j = 0; j < TDim; ++j) {
            rValues[local_index++] = r_node.FastGetSolutionStepValue(r_variables_list[j]->GetTimeDerivative().GetTimeDerivative(), Step);
        }

        rValues[local_index++] = 0.0; // pressure dof
        rValues[local_index++] = r_node.FastGetSolutionStepValue(r_variables_list[TDim+1]->GetTimeDerivative().GetTimeDerivative(), Step);
        rValues[local_index++] = r_node.FastGetSolutionStepValue(r_variables_list[TDim+2]->GetTimeDerivative().GetTimeDerivative(), Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (RansCalculationUtilities::IsWallFunctionActive(this->GetGeometry())) {
        TAdjointConditionData::InitializeCondition(*this, rCurrentProcessInfo);
    }

    this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize and initialize output
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TConditionLocalSize ||
        rLeftHandSideMatrix.size2() != TConditionLocalSize) {
        rLeftHandSideMatrix.resize(TConditionLocalSize, TConditionLocalSize, false);
    }

    rLeftHandSideMatrix.clear();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "TwoEquationTurbulenceModelAdjointCondition::"
                    "CalculateRightHandSide method is not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::CalculateLocalVelocityContribution(
    MatrixType &rDampMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != TConditionLocalSize)
        rRightHandSideVector.resize(TConditionLocalSize, false);

    rRightHandSideVector.clear();

    if (RansCalculationUtilities::IsWallFunctionActive(this->GetGeometry())) {
        auto& r_parent_element = this->GetValue(NEIGHBOUR_ELEMENTS)[0];
        AddFluidResidualsContributions(rRightHandSideVector, r_parent_element, rCurrentProcessInfo);
        AddTurbulenceResidualsContributions(rRightHandSideVector, r_parent_element, rCurrentProcessInfo);
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TConditionLocalSize ||
        rLeftHandSideMatrix.size2() != TConditionLocalSize) {
        rLeftHandSideMatrix.resize(TConditionLocalSize, TConditionLocalSize, false);
    }

    rLeftHandSideMatrix.clear();

    if (RansCalculationUtilities::IsWallFunctionActive(this->GetGeometry())) {
        auto& r_parent_element = this->GetValue(NEIGHBOUR_ELEMENTS)[0];
        AddFluidFirstDerivatives(rLeftHandSideMatrix, r_parent_element, rCurrentProcessInfo);
        AddTurbulenceFirstDerivatives(rLeftHandSideMatrix, r_parent_element, rCurrentProcessInfo);
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TConditionLocalSize ||
        rLeftHandSideMatrix.size2() != TConditionLocalSize) {
        rLeftHandSideMatrix.resize(TConditionLocalSize, TConditionLocalSize, false);
    }

    rLeftHandSideMatrix.clear();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    rMassMatrix.resize(0, 0);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "TwoEquationTurbulenceModelAdjointCondition::"
                    "CalculateDampingMatrix method is not implemented.";

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rSensitivityVariable == SHAPE_SENSITIVITY) {
        if (RansCalculationUtilities::IsWallFunctionActive(this->GetGeometry())) {
            auto& r_parent_elemnet = this->GetValue(NEIGHBOUR_ELEMENTS)[0];
            const IndexType coords_size = TDim * r_parent_elemnet.GetGeometry().PointsNumber();

            if (rOutput.size1() != coords_size || rOutput.size2() != TConditionLocalSize) {
                rOutput.resize(coords_size, TConditionLocalSize, false);
            }

            rOutput.clear();
            BoundedVector<IndexType, TNumNodes> condition_node_parent_index;
            std::vector<IndexType> parent_only_node_indices;
            ComputeParentElementNodesToConditionNodesMap(
                condition_node_parent_index, parent_only_node_indices,
                this->GetGeometry(), r_parent_elemnet.GetGeometry());

            AddFluidShapeDerivatives(rOutput, r_parent_elemnet, condition_node_parent_index, parent_only_node_indices, rCurrentProcessInfo);
            AddTurbulenceShapeDerivatives(rOutput, r_parent_elemnet, condition_node_parent_index, parent_only_node_indices, rCurrentProcessInfo);
        } else {
            constexpr IndexType coords_size = TDim * TNumNodes;

            if (rOutput.size1() != coords_size || rOutput.size2() != TConditionLocalSize) {
                rOutput.resize(coords_size, TConditionLocalSize, false);
            }

            rOutput.clear();
        }
    } else if (rSensitivityVariable == VELOCITY_SENSITIVITY) {
        if (RansCalculationUtilities::IsWallFunctionActive(this->GetGeometry())) {
            auto& r_parent_element = this->GetValue(NEIGHBOUR_ELEMENTS)[0];
            const IndexType coords_size = TDim * r_parent_element.GetGeometry().PointsNumber();

            if (rOutput.size1() != coords_size || rOutput.size2() != TConditionLocalSize) {
                rOutput.resize(coords_size, TConditionLocalSize, false);
            }

            rOutput.clear();

            // compute the lhs which is \partial R \ partial u
            Matrix lhs(TConditionLocalSize, TConditionLocalSize);
            lhs.clear();
            AddFluidFirstDerivatives(lhs, r_parent_element, rCurrentProcessInfo);
            AddTurbulenceFirstDerivatives(lhs, r_parent_element, rCurrentProcessInfo);

            // now need to take a sub set of this matrix
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                    const BoundedVector<double, TConditionLocalSize>& row_values = row(lhs, i_node * TBlockSize + i_dim);
                    AssembleSubVectorToMatrix(rOutput, i_node * TDim + i_dim, 0, row_values);
                }
            }
        } else {
            constexpr IndexType coords_size = TDim * TNumNodes;

            if (rOutput.size1() != coords_size || rOutput.size2() != TConditionLocalSize) {
                rOutput.resize(coords_size, TConditionLocalSize, false);
            }

            rOutput.clear();
        }
    } else {
        KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                     << " not supported." << std::endl;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
std::string TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::Info() const
{
    std::stringstream buffer;
    buffer << "TwoEquationTurbulenceModelAdjointCondition #" << Condition::Id();
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "TwoEquationTurbulenceModelAdjointCondition #" << Condition::Id();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::PrintData(std::ostream& rOStream) const
{
    Condition::pGetGeometry()->PrintData(rOStream);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::AddFluidResidualsContributions(
    Vector& rOutput,
    Element& rParentElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = FluidData::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    this->CalculateGeometryData(Ws, Ns, integration_method);

    using Primal = typename FluidData::Primal;

    typename Primal::Data                   element_data(*this, rParentElement, rCurrentProcessInfo);
    typename Primal::ResidualsContributions residual_contributions(element_data);

    VectorF residual = ZeroVector(TFluidLocalSize);

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const Vector& N = row(Ns, g);
        const double W = Ws[g];

        element_data.CalculateGaussPointData(W, N);
        residual_contributions.AddGaussPointResidualsContributions(residual, W, N);
    }

    AssembleSubVectorToVector(rOutput, 0, residual);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::AddFluidFirstDerivatives(
    MatrixType& rOutput,
    Element& rParentElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = FluidData::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    this->CalculateGeometryData(Ws, Ns, integration_method);

    using Derivatives = typename FluidData::StateDerivatives::FirstDerivatives;

    typename Derivatives::Data                     element_data(*this, rParentElement, rCurrentProcessInfo);
    typename Derivatives::Velocity                 velocity_derivative(element_data);
    typename Derivatives::TurbulenceModelVariable1 turbulence_equation_1_derivative(element_data);
    typename Derivatives::TurbulenceModelVariable2 turbulence_equation_2_derivative(element_data);

    VectorF residual;
    MatrixND dNdXDerivative = ZeroMatrix(TNumNodes, TDim);

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);

        element_data.CalculateGaussPointData(W, N);

        IndexType row = 0;
        for (IndexType c = 0; c < TNumNodes; ++c) {
            for (IndexType k = 0; k < TDim; ++k) {
                velocity_derivative.CalculateGaussPointResidualsDerivativeContributions(residual, c, 999, k, W, N, 0.0, 0.0);
                AssembleSubVectorToMatrix(rOutput, row++, 0, residual);
            }

            // skip pressure derivative
            ++row;

            turbulence_equation_1_derivative.CalculateGaussPointResidualsDerivativeContributions(residual, c, 999, 0, W, N, 0.0, 0.0);
            AssembleSubVectorToMatrix(rOutput, row++, 0, residual);

            turbulence_equation_2_derivative.CalculateGaussPointResidualsDerivativeContributions(residual, c, 999, 0, W, N, 0.0, 0.0);
            AssembleSubVectorToMatrix(rOutput, row++, 0, residual);
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::AddFluidShapeDerivatives(
    Matrix& rOutput,
    Element& rParentElement,
    const BoundedVector<IndexType, TNumNodes>& rConditionNodeParentIndex,
    const std::vector<IndexType>& rParentOnlyNodeIndices,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = FluidData::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    this->CalculateGeometryData(Ws, Ns, integration_method);

    using Derivatives = typename FluidData::SensitivityDerivatives;

    typename Derivatives::Data  element_data(*this, rParentElement, rCurrentProcessInfo);
    typename Derivatives::Shape derivative(element_data);

    VectorF residual_derivative;

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const Vector& N = row(Ns, g);
        const double W = Ws[g];

        element_data.CalculateGaussPointData(W, N);

        IndexType row = 0;

        // calculate derivatives w.r.t. condition nodes
        for (IndexType c = 0; c < TNumNodes; ++c) {
            for (IndexType k = 0; k < TDim; ++k) {
                double detJ_derivative, W_derivative;
                CalculateGeometryDataDerivative(W_derivative, detJ_derivative, c, k, g, integration_method);

                derivative.CalculateGaussPointResidualsDerivativeContributions(residual_derivative, c, rConditionNodeParentIndex[c], k, W, N, W_derivative, detJ_derivative);
                AssembleSubVectorToMatrix(rOutput, row++, 0, residual_derivative);
            }
        }

        // calculate derivatives w.r.t parent element only nodes
        for (IndexType c = 0; c < rParentOnlyNodeIndices.size(); ++c) {
            const IndexType p_index = rParentOnlyNodeIndices[c];
            for (IndexType k = 0; k < TDim; ++k) {
                derivative.CalculateGaussPointResidualsDerivativeContributions(residual_derivative, p_index, k, W, N);
                AssembleSubVectorToMatrix(rOutput, row++, 0, residual_derivative);
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::AddTurbulenceResidualsContributions(
    Vector& rOutput,
    Element& rParentElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& integration_method = TurbulenceModelEquation2Data::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    this->CalculateGeometryData(Ws, Ns, integration_method);

    using Equation1Primal= typename TurbulenceModelEquation2Data::Primal;

    // create data holders for turbulence equations
    typename Equation1Primal::Data eq_2_data(*this, rParentElement, rCurrentProcessInfo);

    // create equation residual data holders
    typename Equation1Primal::ResidualsContributions eq_2_residuals(eq_2_data);

    VectorN residuals = ZeroVector(TNumNodes);

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);

        eq_2_data.CalculateGaussPointData(W, N);
        eq_2_residuals.AddGaussPointResidualsContributions(residuals, W, N);
    }

    AssembleSubVectorToVector(rOutput, TDim + 2, residuals);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::AddTurbulenceFirstDerivatives(
    MatrixType& rOutput,
    Element& rParentElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = TurbulenceModelEquation2Data::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    this->CalculateGeometryData(Ws, Ns, integration_method);

    using Equation2Derivatives = typename TurbulenceModelEquation2Data::StateDerivatives::FirstDerivatives;

    // create data holders for turbulence equations
    typename Equation2Derivatives::Data eq_2_data(*this, rParentElement, rCurrentProcessInfo);

    // create equation derivative data holders
    typename Equation2Derivatives::Velocity                 eq_2_derivative_0(eq_2_data);
    typename Equation2Derivatives::TurbulenceModelVariable1 eq_2_derivative_1(eq_2_data);
    typename Equation2Derivatives::TurbulenceModelVariable2 eq_2_derivative_2(eq_2_data);

    MatrixND dNdX_derivative = ZeroMatrix(TNumNodes, TDim);
    VectorN residual_derivatives;

    IndexType row_index = 0;
    for (IndexType g = 0; g < Ws.size(); ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);

        eq_2_data.CalculateGaussPointData(W, N);

        row_index = 0;
        for (IndexType c = 0; c < TNumNodes; ++c) {

            // add derivatives w.r.t velocity
            for (IndexType k = 0; k < TDim; ++k) {
                eq_2_derivative_0.CalculateGaussPointResidualsDerivativeContributions(residual_derivatives, c, 999, k, W, N, 0, 0);
                AssembleSubVectorToMatrix(rOutput, row_index++, TDim + 2, residual_derivatives);
            }

            // skip derivatives w.r.t. pressure
            ++row_index;

            // add derivatives w.r.t. turbulence variable 1
            eq_2_derivative_1.CalculateGaussPointResidualsDerivativeContributions(residual_derivatives, c, 999, 0, W, N, 0, 0);
            AssembleSubVectorToMatrix(rOutput, row_index++, TDim + 2, residual_derivatives);

            // add derivative w.r.t. turbulence variable 2
            eq_2_derivative_2.CalculateGaussPointResidualsDerivativeContributions(residual_derivatives, c, 999, 0, W, N, 0, 0);
            AssembleSubVectorToMatrix(rOutput, row_index++, TDim + 2, residual_derivatives);
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::AddTurbulenceShapeDerivatives(
    Matrix& rOutput,
    Element& rParentElement,
    const BoundedVector<IndexType, TNumNodes>& rConditionNodeParentIndex,
    const std::vector<IndexType>& rParentOnlyNodeIndices,
    const ProcessInfo& rCurrentProcessInfo)
{
     KRATOS_TRY

    const auto& integration_method = TurbulenceModelEquation2Data::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    this->CalculateGeometryData(Ws, Ns, integration_method);

    using Derivatives = typename TurbulenceModelEquation2Data::SensitivityDerivatives;

    typename Derivatives::Data  element_data(*this, rParentElement, rCurrentProcessInfo);
    typename Derivatives::Shape derivative(element_data);

    VectorN residual_derivative;

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const Vector& N = row(Ns, g);
        const double W = Ws[g];

        element_data.CalculateGaussPointData(W, N);

        IndexType row = 0;

        // calculate derivatives w.r.t. condition nodes
        for (IndexType c = 0; c < TNumNodes; ++c) {
            for (IndexType k = 0; k < TDim; ++k) {
                double detJ_derivative, W_derivative;
                CalculateGeometryDataDerivative(W_derivative, detJ_derivative, c, k, g, integration_method);

                derivative.CalculateGaussPointResidualsDerivativeContributions(residual_derivative, c, rConditionNodeParentIndex[c], k, W, N, W_derivative, detJ_derivative);
                AssembleSubVectorToMatrix(rOutput, row++, TDim + 2, residual_derivative);
            }
        }

        // calculate derivatives w.r.t parent element only nodes
        for (IndexType c = 0; c < rParentOnlyNodeIndices.size(); ++c) {
            const IndexType p_index = rParentOnlyNodeIndices[c];
            for (IndexType k = 0; k < TDim; ++k) {
                derivative.CalculateGaussPointResidualsDerivativeContributions(residual_derivative, p_index, k, W, N);
                AssembleSubVectorToMatrix(rOutput, row++, TDim + 2, residual_derivative);
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::ComputeParentElementNodesToConditionNodesMap(
    BoundedVector<IndexType, TNumNodes>& rConditionNodeParentIndex,
    std::vector<IndexType>& rParentOnlyNodeIndices,
    const GeometryType& rConditionGeometry,
    const GeometryType& rParentElementGeometry) const
{
    const IndexType parent_element_number_of_nodes = rParentElementGeometry.PointsNumber();

    rParentOnlyNodeIndices.resize(parent_element_number_of_nodes - TNumNodes);

    // collect condition and parent element common node information
    IndexType parent_nodes = 0;
    for (IndexType p_i = 0; p_i < parent_element_number_of_nodes; ++p_i) {
        IndexType c_i;
        for (c_i = 0; c_i < TNumNodes; ++c_i) {
            if (rConditionGeometry[c_i].Id() == rParentElementGeometry[p_i].Id()) {
                rConditionNodeParentIndex[c_i] = p_i;
                break;
            }
        }

        // if the parent element node is not found in the condition
        if (c_i == TNumNodes) {
            rParentOnlyNodeIndices[parent_nodes++] = p_i;
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::CalculateGeometryData(
    Vector& rGaussWeights,
    Matrix& rNContainer,
    const GeometryData::IntegrationMethod& rIntegrationMethod) const
{
    RansCalculationUtilities::CalculateConditionGeometryData(this->GetGeometry(), rIntegrationMethod, rGaussWeights, rNContainer);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointConditionData>
void TwoEquationTurbulenceModelAdjointCondition<TDim, TNumNodes, TAdjointConditionData>::CalculateGeometryDataDerivative(
    double& WDerivative,
    double& DetJDerivative,
    const IndexType NodeIndex,
    const IndexType DirectionIndex,
    const IndexType GaussPointIndex,
    const GeometryData::IntegrationMethod& rIntegrationMethod) const
{
    const auto& integration_points = this->GetGeometry().IntegrationPoints(rIntegrationMethod);

    const double domain_size_derivative = RansAdjointUtilities::GeometricalDerivatives<TDim, TNumNodes>::DomainSizeDerivative(this->GetGeometry(), NodeIndex, DirectionIndex);

    // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
    DetJDerivative = (TDim == 2) ? 0.5 * domain_size_derivative : 2.0 * domain_size_derivative;
    WDerivative = DetJDerivative * integration_points[GaussPointIndex].Weight();
}

// template instantiations
// k-epsilon
template class TwoEquationTurbulenceModelAdjointCondition<2, 2, KEpsilonWallConditionData::VMSMonolithicKBasedEpsilonKBasedWallConditionData<2, 2>>;
template class TwoEquationTurbulenceModelAdjointCondition<3, 3, KEpsilonWallConditionData::VMSMonolithicKBasedEpsilonKBasedWallConditionData<3, 3>>;

template class TwoEquationTurbulenceModelAdjointCondition<2, 2, KEpsilonWallConditionData::VMSMonolithicKBasedEpsilonUBasedWallConditionData<2, 2>>;
template class TwoEquationTurbulenceModelAdjointCondition<3, 3, KEpsilonWallConditionData::VMSMonolithicKBasedEpsilonUBasedWallConditionData<3, 3>>;

// k-omega
template class TwoEquationTurbulenceModelAdjointCondition<2, 2, KOmegaWallConditionData::VMSMonolithicKBasedOmegaKBasedWallConditionData<2, 2>>;
template class TwoEquationTurbulenceModelAdjointCondition<3, 3, KOmegaWallConditionData::VMSMonolithicKBasedOmegaKBasedWallConditionData<3, 3>>;

template class TwoEquationTurbulenceModelAdjointCondition<2, 2, KOmegaWallConditionData::VMSMonolithicKBasedOmegaUBasedWallConditionData<2, 2>>;
template class TwoEquationTurbulenceModelAdjointCondition<3, 3, KOmegaWallConditionData::VMSMonolithicKBasedOmegaUBasedWallConditionData<3, 3>>;

// k-omega-sst
template class TwoEquationTurbulenceModelAdjointCondition<2, 2, KOmegaSSTWallConditionData::VMSMonolithicKBasedOmegaKBasedWallConditionData<2, 2>>;
template class TwoEquationTurbulenceModelAdjointCondition<3, 3, KOmegaSSTWallConditionData::VMSMonolithicKBasedOmegaKBasedWallConditionData<3, 3>>;

template class TwoEquationTurbulenceModelAdjointCondition<2, 2, KOmegaSSTWallConditionData::VMSMonolithicKBasedOmegaUBasedWallConditionData<2, 2>>;
template class TwoEquationTurbulenceModelAdjointCondition<3, 3, KOmegaSSTWallConditionData::VMSMonolithicKBasedOmegaUBasedWallConditionData<3, 3>>;

} // namespace Kratos