//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Kratos Team
//

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "future/includes/adjoint_derivatives_interface.h"
#include "includes/node.h"

namespace Kratos::Future {

///@addtogroup KratosCore
///@{

/**
 * @brief Opt-in implementation of AdjointDerivativesInterface for linear entities (elements or conditions).
 * @details Usage: class MyElement : public Element, public LinearAdjointDerivatives<MyElement>
 *          Derivatives are taken from the existing CalculateLeftHandSide, CalculateDampingMatrix and
 *          CalculateMassMatrix, so it is only valid when those are the exact linear matrices (no Picard or secant ones).
 *          The adjoint DOFs are the adjoint of the primal DOFs, in the order of the primal GetDofList.
 *          Only historical nodal data gives non-zero derivatives, with respect to the DOF variable (Stiffness),
 *          its first time derivative (Damping) and its second time derivative (Inertia).
 * @tparam TEntity Entity type, which must provide GetDofList, GetGeometry and the matrix calculation methods
 */
template<class TEntity>
class LinearAdjointDerivatives : public AdjointDerivativesInterface
{
public:
    ///@name Operations
    ///@{

    void GetAdjointDofList(
        const std::vector<DerivativeVariable>& rVariables,
        DofsVectorType& rAdjointDofs,
        const ProcessInfo& rProcessInfo) const override
    {
        KRATOS_TRY

        DofsVectorType primal_dofs;
        GetEntity().GetDofList(primal_dofs, rProcessInfo);

        rAdjointDofs.clear();
        rAdjointDofs.reserve(primal_dofs.size());
        for (const auto& p_dof : primal_dofs) {
            const auto& r_variable = DofVariable(*p_dof);
            KRATOS_ERROR_IF_NOT(r_variable.HasAdjointVariable())
                << "Primal DOF variable " << r_variable.Name() << " has no adjoint variable." << std::endl;
            const auto& r_node = FindNode(p_dof->Id());
            KRATOS_ERROR_IF_NOT(r_node.HasDofFor(r_variable.GetAdjointVariable()))
                << "Node " << r_node.Id() << " has no DOF for " << r_variable.GetAdjointVariable().Name() << std::endl;
            rAdjointDofs.push_back(r_node.pGetDof(r_variable.GetAdjointVariable()));
        }

        KRATOS_CATCH("")
    }

    void GetAdjointEquationIdVector(
        const std::vector<DerivativeVariable>& rVariables,
        EquationIdVectorType& rEquationIds,
        const ProcessInfo& rProcessInfo) const override
    {
        KRATOS_TRY

        DofsVectorType adjoint_dofs;
        GetAdjointDofList(rVariables, adjoint_dofs, rProcessInfo);
        rEquationIds.resize(adjoint_dofs.size());
        for (std::size_t i = 0; i < adjoint_dofs.size(); ++i) {
            rEquationIds[i] = adjoint_dofs[i]->EquationId();
        }

        KRATOS_CATCH("")
    }

    void GetAdjointValues(
        const std::vector<DerivativeVariable>& rVariables,
        Vector& rValues,
        const ProcessInfo& rProcessInfo,
        const std::size_t StepIndex = 0) const override
    {
        GetValues(0, rVariables, rValues, rProcessInfo, StepIndex);
    }

    void GetAdjointFirstDerivativeValues(
        const std::vector<DerivativeVariable>& rVariables,
        Vector& rValues,
        const ProcessInfo& rProcessInfo,
        const std::size_t StepIndex = 0) const override
    {
        GetValues(1, rVariables, rValues, rProcessInfo, StepIndex);
    }

    void GetAdjointSecondDerivativeValues(
        const std::vector<DerivativeVariable>& rVariables,
        Vector& rValues,
        const ProcessInfo& rProcessInfo,
        const std::size_t StepIndex = 0) const override
    {
        GetValues(2, rVariables, rValues, rProcessInfo, StepIndex);
    }

    void CalculateResidualTermDerivatives(
        const TermDerivativesType& rTermDerivatives,
        const std::vector<DerivativeVariable>& rVariables,
        const ProcessInfo& rProcessInfo,
        const std::size_t StepIndex = 0) const override // The linear matrices are the same for every step
    {
        KRATOS_TRY

        DofsVectorType primal_dofs;
        GetEntity().GetDofList(primal_dofs, rProcessInfo);
        const std::size_t number_of_nodes = GetEntity().GetGeometry().size();

        // The linear matrices are calculated only once and only if a term needs them
        TEntity& r_entity = const_cast<TEntity&>(GetEntity());
        Matrix matrices[3];
        bool computed[3] = {false, false, false};

        std::vector<std::size_t> offsets(rVariables.size() + 1, 0);
        for (std::size_t v = 0; v < rVariables.size(); ++v) {
            offsets[v + 1] = offsets[v] + (IsNodal(rVariables[v]) ? number_of_nodes : 1);
        }

        for (const auto& [term, r_matrix_wrapper] : rTermDerivatives) {
            Matrix& r_output = r_matrix_wrapper.get();
            r_output.resize(0, 0, false);
            const int order = TimeOrder(term);
            if (order < 0) {
                continue; // Force and Constraint have no derivative from the linear matrices
            }
            if (!computed[order]) {
                if (order == 0) {
                    r_entity.CalculateLeftHandSide(matrices[0], rProcessInfo);
                } else if (order == 1) {
                    r_entity.CalculateDampingMatrix(matrices[1], rProcessInfo);
                } else {
                    r_entity.CalculateMassMatrix(matrices[2], rProcessInfo);
                }
                computed[order] = true;
            }
            const Matrix& r_matrix = matrices[order];

            r_output = ZeroMatrix(primal_dofs.size(), offsets.back());
            for (std::size_t v = 0; v < rVariables.size(); ++v) {
                const auto& r_requested = rVariables[v];
                if (r_requested.GetLocation() != Globals::DataLocation::NodeHistorical) {
                    continue;
                }
                for (std::size_t j = 0; j < primal_dofs.size(); ++j) {
                    const auto* p_state = StateVariable(DofVariable(*primal_dofs[j]), order);
                    if (p_state == nullptr || !Matches(r_requested, *p_state)) {
                        continue;
                    }
                    const std::size_t column = offsets[v] + NodeIndex(primal_dofs[j]->Id());
                    for (std::size_t i = 0; i < primal_dofs.size(); ++i) {
                        r_output(i, column) = r_matrix(i, j);
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    void GetValues(
        const int Order,
        const std::vector<DerivativeVariable>& rVariables,
        Vector& rValues,
        const ProcessInfo& rProcessInfo,
        const std::size_t StepIndex) const
    {
        KRATOS_TRY

        DofsVectorType adjoint_dofs;
        GetAdjointDofList(rVariables, adjoint_dofs, rProcessInfo);
        rValues.resize(adjoint_dofs.size(), false);
        for (std::size_t i = 0; i < adjoint_dofs.size(); ++i) {
            const auto* p_variable = StateVariable(DofVariable(*adjoint_dofs[i]), Order);
            rValues[i] = p_variable == nullptr ? 0.0 : FindNode(adjoint_dofs[i]->Id()).GetSolutionStepValue(*p_variable, StepIndex);
        }

        KRATOS_CATCH("")
    }

    const TEntity& GetEntity() const
    {
        return static_cast<const TEntity&>(*this);
    }

    // 0: Stiffness, 1: Damping, 2: Inertia, -1: no linear derivative
    static int TimeOrder(const ResidualTerm Term)
    {
        switch (Term) {
            case ResidualTerm::Stiffness: return 0;
            case ResidualTerm::Damping: return 1;
            case ResidualTerm::Inertia: return 2;
            default: return -1;
        }
    }

    static const Variable<double>& DofVariable(const Dof<double>& rDof)
    {
        return static_cast<const Variable<double>&>(rDof.GetVariable());
    }

    // Variable multiplied by the matrix of the term, or nullptr if it does not exist
    static const Variable<double>* StateVariable(const Variable<double>& rDofVariable, const int Order)
    {
        const Variable<double>* p_variable = &rDofVariable;
        for (int i = 0; i < Order; ++i) {
            if (!p_variable->HasTimeDerivative()) {
                return nullptr;
            }
            p_variable = &p_variable->GetTimeDerivative();
        }
        return p_variable;
    }

    static bool IsNodal(const DerivativeVariable& rRequested)
    {
        return rRequested.GetLocation() == Globals::DataLocation::NodeHistorical
            || rRequested.GetLocation() == Globals::DataLocation::NodeNonHistorical;
    }

    // Checks if the requested variable refers to the scalar (or component) variable
    static bool Matches(const DerivativeVariable& rRequested, const Variable<double>& rActual)
    {
        return rRequested.GetVariable().Key() == rActual.Key() && rRequested.GetComponentIndex() == 0;
    }

    std::size_t NodeIndex(const std::size_t NodeId) const
    {
        const auto& r_geometry = GetEntity().GetGeometry();
        for (std::size_t i = 0; i < r_geometry.size(); ++i) {
            if (r_geometry[i].Id() == NodeId) {
                return i;
            }
        }
        KRATOS_ERROR << "Node " << NodeId << " not found in the geometry of the entity." << std::endl;
    }

    const Node& FindNode(const std::size_t NodeId) const
    {
        return GetEntity().GetGeometry()[NodeIndex(NodeId)];
    }

    ///@}
};

///@}

} // namespace Kratos::Future
