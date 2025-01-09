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

#pragma once

// System includes
#include <atomic>
#include <vector>
#include <tuple>

// External includes

// Project includes
#include "containers/variable.h"
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/process_info.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"

namespace Kratos {

/**
 * @brief Base EntityAdjointExtension class
 * @details This class is used to extend the primal entity with information
 *          on how to compute partial derivatives of the primal terms in the residual for set of given variables.
 *
 * @tparam TEntity  Type of the entity (i.e. Element/Condition/...)
 */
template <class TEntity>
class EntityAdjointExtension {
public:
    ///@name Type definitions
    ///@{

    enum Term
    {
        MASS,
        DAMPING,
        STIFFNESS,
        FORCE
    };

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EntityAdjointExtension);

    using DofType = Dof<double>;

    using EquationIdVectorType = std::vector<std::size_t>;

    using DofsVectorType = std::vector<DofType::Pointer>;

    using SensitivityMatrixType = std::vector<Matrix>;

    using TermSensitivityMatrixType = std::pair<Term, SensitivityMatrixType*>;

    ///@}
    ///@name Life cycle
    ///@{

    EntityAdjointExtension() = default;

    virtual ~EntityAdjointExtension() = default;

    ///@}
    ///@name Adjoint variable operations
    ///@{

    virtual void AdjointEquationIdVector(
        EquationIdVectorType& rEquationIds,
        const TEntity& rEntity,
        const ProcessInfo& rProcessInfo) const {}

    virtual void GetAdjointDofList(
        DofsVectorType& rDofList,
        const TEntity& rEntity,
        const ProcessInfo& rProcessInfo) const {}

    virtual void GetAdjointValuesVector(
        Vector& rValues,
        const TEntity& rEntity,
        const int Step = 0) const  {}

    virtual void GetAdjointFirstDerivativesVector(
        Vector& rValues,
        const TEntity& rEntity,
        const int Step = 0) const  {}

    virtual void GetAdjointSecondDerivativesVector(
        Vector& rValues,
        const TEntity& rEntity,
        const int Step = 0) const  {}

    virtual void SetAdjointFirstDerivativesVector(
        TEntity& rEntity,
        const Vector& rValues,
        const int Step = 0) const {}

    virtual void SetAdjointSecondDerivativesVector(
        TEntity& rEntity,
        const Vector& rValues,
        const int Step = 0) const {}

    ///@}
    ///@name State variable partial derivative operations
    ///@{

    /**
     * @brief Adds a state variable for which the partial derivatives are required.
     * @details The order of added state variables decides the order of the partial derivatives given in the
     *          third order tensor.
     * @param rVariable The double state variable to be added.
     */
    virtual void AddStateVariable(const Variable<double>& rVariable) {}

    /**
     * @brief Calculates the partial derivatives w.r.t. given list of state variables.
     * @details This method is used to compute the partial derivatives on multiple terms (i.e. MASS, DAMPING)
     *          w.r.t. given list of state variables. A list of terms is provided to take advantage of computing
     *          constants and common things once for all the terms if required. The implementation is of this
     *          method is not restricted, hence developers can implement the method how ever they please.
     *
     *          the output of SensitivityMatrixType will have the following format for the mass matrix
     *          where \f(m_{ij}\f) is the \f(i^{th}\f) row and \f(j^{th}\f) column value and computes
     *          the derivative w.r.t. state variable \f(\phi_k\f) then:
     *          \f[
     *              Output(k, i, j) = \frac{\partial m_{ij}}{\partial \phi_k}
     *          \f]
     *
     * @param rOutput           List of pair containing the term and the Output(k, i, j) tensor for each term.
     * @param rEntity           Entity to be used to compute the derivatives of the residual.
     * @param rProcessInfo      Current Processinfo.
     */
    virtual void CalculatePartialStateVariableDerivatives(
        std::vector<TermSensitivityMatrixType>&& rOutput,
        TEntity& rEntity,
        const ProcessInfo& rProcessInfo) const {}

    /**
     * @brief Adds a design variable for which the partial derivatives are required.
     * @details The order of added design variables decides the order of the partial derivatives given in the
     *          third order tensor.
     * @param rVariable The double design variable to be added.
     */
    virtual void AddDesignVariable(const Variable<double>& rVariable) {}

    /**
     * @brief Calculates the partial derivatives w.r.t. given list of design variables.
     * @details This method is used to compute the partial derivatives on multiple terms (i.e. MASS, DAMPING)
     *          w.r.t. given list of design variables. A list of terms is provided to take advantage of computing
     *          constants and common things once for all the terms if required. The implementation is of this
     *          method is not restricted, hence developers can implement the method how ever they please.
     *
     *          the output of SensitivityMatrixType will have the following format for the mass matrix
     *          where \f(m_{ij}\f) is the \f(i^{th}\f) row and \f(j^{th}\f) column value and computes
     *          the derivative w.r.t. design variable \f(\phi_k\f) then:
     *          \f[
     *              Output(k, i, j) = \frac{\partial m_{ij}}{\partial \phi_k}
     *          \f]
     *
     * @param rOutput           List of pair containing the term and the Output(k, i, j) tensor for each term.
     * @param rEntity           Entity to be used to compute the derivatives of the residual.
     * @param rProcessInfo      Current Processinfo.
     */
    virtual void CalculatePartialDesignVariableDerivatives(
        std::vector<TermSensitivityMatrixType>&& rOutput,
        TEntity& rEntity,
        const ProcessInfo& rProcessInfo) const {}

    ///@}
    ///@name Inquiry
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "EntityAdjointExtension";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}

private:
    ///@name Private Operators
    ///@{

    //*********************************************
    // this block is needed for refcounting
    mutable std::atomic<int> mReferenceCounter{0};

    friend void intrusive_ptr_add_ref(const EntityAdjointExtension<TEntity>* x)
    {
        x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
    }

    friend void intrusive_ptr_release(const EntityAdjointExtension<TEntity>* x)
    {
        if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
            std::atomic_thread_fence(std::memory_order_acquire);
            delete x;
        }
    }
    //*********************************************

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }

    ///@}
};

///@name Input and output
///@{

/// output stream function
template<class TEntity>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const EntityAdjointExtension<TEntity>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}
///@}

} // namespace Kratos