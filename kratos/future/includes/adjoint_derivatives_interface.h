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
#include <functional>
#include <utility>
#include <vector>

// External includes

// Project includes
#include "containers/variable_data.h"
#include "containers/variable.h"
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/global_variables.h"
#include "includes/ublas_interface.h"
#include "includes/process_info.h"

namespace Kratos::Future {

///@addtogroup KratosCore
///@{

/**
 * @brief Terms of the residual R = F - M*a - C*v - K*u (plus constraints) whose derivatives can be requested.
 */
enum class ResidualTerm
{
    Force,      /// F
    Inertia,    /// M*a
    Damping,    /// C*v
    Stiffness,  /// K*u
    Constraint  /// g = s - T*m - c (master-slave and multi-freedom constraints)
};

/**
 * @brief Describes the variable (component) with respect to which derivatives are requested.
 * @details The variable is always a primal variable. Adjoint variables are found through Variable::GetAdjointVariable.
 *          Components of array variables are given as component variables (e.g. DISPLACEMENT_X).
 */
class DerivativeVariable
{
public:
    /// Scalar variable or component variable
    DerivativeVariable(const Variable<double>& rVariable, const Globals::DataLocation Location = Globals::DataLocation::NodeHistorical)
        : mpVariable(&rVariable), mLocation(Location)
    {
    }

    /// Entry of a Vector variable
    DerivativeVariable(const Variable<Vector>& rVariable, const std::size_t ComponentIndex, const Globals::DataLocation Location = Globals::DataLocation::NodeHistorical)
        : mpVariable(&rVariable), mLocation(Location), mComponentIndex(ComponentIndex)
    {
    }

    /// Entry of a Matrix variable (flat, row-major index)
    DerivativeVariable(const Variable<Matrix>& rVariable, const std::size_t ComponentIndex, const Globals::DataLocation Location = Globals::DataLocation::NodeHistorical)
        : mpVariable(&rVariable), mLocation(Location), mComponentIndex(ComponentIndex)
    {
    }

    const VariableData& GetVariable() const noexcept { return *mpVariable; }

    Globals::DataLocation GetLocation() const noexcept { return mLocation; }

    /// Flat index into the Vector or Matrix variable, 0 for scalar and component variables
    std::size_t GetComponentIndex() const noexcept { return mComponentIndex; }

private:
    const VariableData* mpVariable;
    const Globals::DataLocation mLocation;
    const std::size_t mComponentIndex = 0;
};

/**
 * @brief Interface to be inherited by elements, conditions, master-slave constraints and multi-freedom constraints
 *        to provide the derivatives needed to build and solve the adjoint problem.
 * @details There are no default implementations. Every method raises an error unless overridden.
 *          The order of the adjoint DOFs is defined by the implementor. GetAdjointDofList and
 *          GetAdjointEquationIdVector must return the same entries in the same order, and that order
 *          is the one used by the rows of the matrices of CalculateResidualTermDerivatives.
 */
class KRATOS_API(KRATOS_CORE) AdjointDerivativesInterface
{
public:
    ///@name Type Definitions
    ///@{

    using DofType = Dof<double>;
    using DofsVectorType = std::vector<DofType::Pointer>;
    using EquationIdVectorType = std::vector<std::size_t>;

    /// Residual terms paired with the matrices where their derivatives are stored. Example: {{ResidualTerm::Inertia, mass}, {ResidualTerm::Damping, damping}}
    using TermDerivativesType = std::vector<std::pair<ResidualTerm, std::reference_wrapper<Matrix>>>;

    ///@}
    ///@name Life Cycle
    ///@{

    virtual ~AdjointDerivativesInterface() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Gets the adjoint DOFs of this entity for the given set of derivative variables.
     * @details The list is not necessarily the adjoint of the primal DOFs. In coupled problems it can contain
     *          adjoint DOFs of fields not owned by this entity. The order is defined by the implementor.
     */
    virtual void GetAdjointDofList(
        const std::vector<DerivativeVariable>& rVariables,
        DofsVectorType& rAdjointDofs,
        const ProcessInfo& rProcessInfo) const
    {
        KRATOS_ERROR << "GetAdjointDofList is not implemented." << std::endl;
    }

    /**
     * @brief Gets the equation ids of the adjoint DOFs, in the same order as GetAdjointDofList.
     */
    virtual void GetAdjointEquationIdVector(
        const std::vector<DerivativeVariable>& rVariables,
        EquationIdVectorType& rEquationIds,
        const ProcessInfo& rProcessInfo) const
    {
        KRATOS_ERROR << "GetAdjointEquationIdVector is not implemented." << std::endl;
    }

    /**
     * @brief Gets the adjoint values of the adjoint DOFs, read from the nodal data of the step StepIndex steps back (0 is the current step).
     * @details One entry per adjoint DOF, in the order of GetAdjointDofList. Entries without a time derivative are 0.
     */
    virtual void GetAdjointValues(
        const std::vector<DerivativeVariable>& rVariables,
        Vector& rValues,
        const ProcessInfo& rProcessInfo,
        const std::size_t StepIndex = 0) const
    {
        KRATOS_ERROR << "GetAdjointValues is not implemented." << std::endl;
    }

    /**
     * @brief Gets the first time derivatives of the adjoint values of the adjoint DOFs, read from the nodal data of the step StepIndex steps back (0 is the current step).
     * @details One entry per adjoint DOF, in the order of GetAdjointDofList. Entries without a time derivative are 0.
     */
    virtual void GetAdjointFirstDerivativeValues(
        const std::vector<DerivativeVariable>& rVariables,
        Vector& rValues,
        const ProcessInfo& rProcessInfo,
        const std::size_t StepIndex = 0) const
    {
        KRATOS_ERROR << "GetAdjointFirstDerivativeValues is not implemented." << std::endl;
    }

    /**
     * @brief Gets the second time derivatives of the adjoint values of the adjoint DOFs, read from the nodal data of the step StepIndex steps back (0 is the current step).
     * @details One entry per adjoint DOF, in the order of GetAdjointDofList. Entries without a time derivative are 0.
     */
    virtual void GetAdjointSecondDerivativeValues(
        const std::vector<DerivativeVariable>& rVariables,
        Vector& rValues,
        const ProcessInfo& rProcessInfo,
        const std::size_t StepIndex = 0) const
    {
        KRATOS_ERROR << "GetAdjointSecondDerivativeValues is not implemented." << std::endl;
    }

    /**
     * @brief Calculates the derivatives of the given residual terms with respect to the given variables.
     * @details For each (term, matrix) pair in rTermDerivatives the matrix is filled with:
     *          rows    = adjoint DOFs of this entity (order of GetAdjointDofList),
     *          columns = the entries of all rVariables on this entity, concatenated in the order of rVariables.
     *                    Each variable occupies 1 column for each node for nodal data (historical or not), and 1 column otherwise.
     *          An empty (0x0) matrix means identically zero.
     *          StepIndex selects the residual: 0 is the residual of the current step, k is the residual of k steps back.
     *          The derivatives are taken with respect to rVariables of that same step. The coupling between steps
     *          comes from the time integration and is handled by the adjoint scheme. The primal state of that step
     *          is expected to be available in the entity (e.g. restored by restart).
     */
    virtual void CalculateResidualTermDerivatives(
        const TermDerivativesType& rTermDerivatives,
        const std::vector<DerivativeVariable>& rVariables,
        const ProcessInfo& rProcessInfo,
        const std::size_t StepIndex = 0) const
    {
        KRATOS_ERROR << "CalculateResidualTermDerivatives is not implemented." << std::endl;
    }

    ///@}
};

///@}

} // namespace Kratos::Future
