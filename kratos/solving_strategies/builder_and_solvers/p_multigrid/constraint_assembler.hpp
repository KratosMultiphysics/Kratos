//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/status_stream.hpp" // PMGStatusStream
#include "containers/data_value_container.h" // DataValueContainer
#include "containers/pointer_vector_set.h" // PointerVectorSet
#include "includes/master_slave_constraint.h" // MasterSlaveConstraint
#include "includes/indexed_object.h" // IndexedObject
#include "includes/dof.h" // Dof
#include "includes/variables.h" // IDENTIFIER
#include "includes/smart_pointers.h" // KRATOS_CLASS_POINTER_DEFINITION


namespace Kratos {


/// @brief Enum class representing imposition methods for multifreedom constraints.
enum class ConstraintImposition
{
    None                = 0,
    MasterSlave         = 1,
    Lagrange            = 2,
    AugmentedLagrange   = 3,
    Penalty             = 4
}; // enum class ConstraintImposition


/// @brief Interface for assembling and imposing multifreedom constraints.
template <class TSparse, class TDense>
class ConstraintAssembler : public DataValueContainer
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ConstraintAssembler);

    using DofSet = PointerVectorSet<Dof<typename TDense::DataType>>;

    using ConstraintArray = PointerVectorSet<MasterSlaveConstraint, IndexedObject>;

    ConstraintAssembler() noexcept
        : ConstraintAssembler(ConstraintImposition::MasterSlave, "unnamed")
    {}

    ConstraintAssembler(ConstraintImposition Method,
                        std::string&& rInstanceName)
        : DataValueContainer(),
          mRelationMatrix(),
          mConstraintGapVector(),
          mName(std::move(rInstanceName))
    {
        std::string method_name;

        switch (Method) {
            case ConstraintImposition::None:
                method_name = "none";
                break;
            case ConstraintImposition::MasterSlave:
                method_name = "master_slave";
                break;
            case ConstraintImposition::AugmentedLagrange:
                method_name = "augmented_lagrange";
                break;
            default:
                // Other imposition methods are not implemented yet.
                KRATOS_ERROR << "Unsupported constraint imposition: " << (int)Method;
        } // switch Method

        this->SetValue(ConstraintAssembler::GetImpositionVariable(), method_name);
    }

    ConstraintAssembler(ConstraintAssembler&&) noexcept = default;

    /// @details Define an overriding virtual destructor to ensure compile time errors
    //           if the base class' destructor turns non-virtual.
    virtual ~ConstraintAssembler() override = default;

    ConstraintAssembler& operator=(ConstraintAssembler&&) noexcept = default;

    /// @brief Allocate memory for the constraint gap vector and relation matrix, and compute its topology.
    /// @details This function is responsible for large memory allocations, as well as computing
    ///          the sparsity pattern of the relation matrix. It must also modify the provided
    ///          left hand side matrix, solution vector, right hand side vector, and DoF list such that
    ///          these containers will not require reallocation during later stages of the solution process.
    /// @param rConstraints Constraint set of the related @ref ModelPart.
    /// @param rProcessInfo Current @ref ProcessInfo of the related @ref ModelPart.
    /// @param rLhs Unconstrained left hand side matrix' topology.
    /// @param rDofSet Unconstrained set of @ref Dof "DoFs".
    /// @note This function should be invoked @b after the unconstrained system is allocated, but @b before
    ///       it is assembled.
    virtual void Allocate(const ConstraintArray& rConstraints,
                          const ProcessInfo& rProcessInfo,
                          typename TSparse::MatrixType& rLhs,
                          typename TSparse::VectorType& rSolution,
                          typename TSparse::VectorType& rRhs,
                          DofSet& rDofSet)
    {
    }

    /// @brief Compute and assemble constraint contributions into the preallocated relation matrix and constraint gap vector.
    /// @details This function is responsible for computing the entries of the relation matrix
    ///          as well as the constraint gap vector.
    /// @param rConstraints Constraint set of the related @ref ModelPart.
    /// @param rProcessInfo @ref ProcessInfo of the related @ref ModelPart.
    /// @param rDofSet Unconstrained set of @ref Dof "DoFs".
    /// @param AssembleLhs Indicates whether to assemble data structures necessary for updating
    ///                    the left hand side matrix of the unconstrained system.
    /// @param AssembleRhs Indicates whether to assemble data structures necessary for updating
    ///                    the right hand side vector of the unconstrained system.
    /// @note This function must be preceded by a call to @ref ConstraintAssembler::Allocate, and should not make large scale
    ///       reallocations.
    virtual void Assemble(const ConstraintArray& rConstraints,
                          const ProcessInfo& rProcessInfo,
                          DofSet& rDofSet,
                          const bool AssembleLhs,
                          const bool AssembleRhs)
    {
    }

    /// @brief Prepare the linear system for the solution loop.
    /// @details This function is supposed to perform tasks on the linear system
    ///          that are required only once, before calls to the linear solver begin.
    ///          Constraint imposition methods that do not require a solution loop
    ///          (for example, master-slave elimination one-shots the constraints)
    ///          should manipulate the system here. If the set of @ref Dof "DoFs" has
    ///          to be changed, it should also be carried out here.
    /// @param rLhs Unconstrained left hand side matrix with topology to accomodate constraint imposition.
    /// @param rRhs Unconstrained right hand side vector with space to accomodate constraint imposition.
    /// @param itDofBegin Iterator pointing to the first @ref Dof "DoF" of the unconstrained system.
    /// @param itDofEnd Sentinel of the unconstrained system's array of @ref Dof "DoFs".
    virtual void Initialize(typename TSparse::MatrixType& rLhs,
                            typename TSparse::VectorType& rRhs,
                            typename DofSet::iterator itDofBegin,
                            typename DofSet::iterator itDofEnd)
    {
    }

    /// @brief Manipulate the linear system before invoking the linear solver in the solution loop's current iteration.
    /// @param rLhs Left hand side matrix.
    /// @param rSolution Unconverged solution vector.
    /// @param rRhs Right hand side vector.
    virtual void InitializeSolutionStep(typename TSparse::MatrixType& rLhs,
                                        typename TSparse::VectorType& rSolution,
                                        typename TSparse::VectorType& rRhs)
    {
    }

    /// @brief Perform constraint-related tasks after invoking the linear solver in the current iteration of the solution loop.
    /// @details This function is supposed to evaluate the convergence of constraint imposition,
    ///          decide whether to request more iterations in the solution loop.
    /// @param rLhs Constrained left hand side matrix.
    /// @param rSolution Converged solution vector.
    /// @param rRhs Constrained right hand side vector.
    /// @param rReport Status information on the solution loop.
    /// @warning The solution loop will continue indefinitely unless this function eventually
    ///          returns @p true.
    virtual bool FinalizeSolutionStep(typename TSparse::MatrixType& rLhs,
                                      typename TSparse::VectorType& rSolution,
                                      typename TSparse::VectorType& rRhs,
                                      PMGStatusStream::Report& rReport) = 0;

    /// @brief Perform tasks related to constraint imposition after constraints converged.
    /// @param rLhs Constrained left hand side matrix.
    /// @param rSolution Converged solution vector.
    /// @param rRhs Constrained right hand side vector.
    /// @param rDofSet Constrained set of @ref Dof "DoFs".
    virtual void Finalize(typename TSparse::MatrixType& rLhs,
                          typename TSparse::VectorType& rSolution,
                          typename TSparse::VectorType& rRhs,
                          DofSet& rDofSet)
    {
    }

    /// @brief Release memory related to the linear system and constraints.
    /// @details Derived classes must call the @ref ConstraintAssembler::Clear "Clear"
    ///          function of their parents at some point.
    virtual void Clear()
    {
        mRelationMatrix = typename TSparse::MatrixType();
        mConstraintGapVector = typename TSparse::VectorType();
    }

    ConstraintImposition GetImposition() const
    {
        const std::string method_name = this->GetValue(ConstraintAssembler::GetImpositionVariable());
        if (method_name == "none") {
            return ConstraintImposition::None;
        } else if (method_name == "master_slave") {
            return ConstraintImposition::MasterSlave;
        } else if (method_name == "augmented_lagrange") {
            return ConstraintImposition::AugmentedLagrange;
        } else {
            // Other imposition methods are not implemented yet.
            KRATOS_ERROR << "Unsupported constraint imposition: \"" << method_name << "\"";
        }
    }

    const typename TSparse::MatrixType& GetRelationMatrix() const noexcept
    {
        return mRelationMatrix;
    }

    const typename TSparse::VectorType& GetConstraintGapVector() const noexcept
    {
        return mConstraintGapVector;
    }

    const std::string& GetName() const noexcept
    {
        return mName;
    }

protected:
    // PGrid inherits an assembled unconstrained system, as well as an assembled relation matrix,
    // and constructs a restriction operator it can use to directly compute the system and relation
    // matrix at its own level. As a result, it does not need to do allocation and assembly, which
    // means it should be able to directly define the relation matrix and constraint gap vector
    // of its own constraint assembler.
    // => long story short, PGrid needs access to the protected members of this class.
    template <class S, class D>
    friend class PGrid;

    typename TSparse::MatrixType& GetRelationMatrix() noexcept
    {
        return mRelationMatrix;
    }

    typename TSparse::VectorType& GetConstraintGapVector() noexcept
    {
        return mConstraintGapVector;
    }

private:
    ConstraintAssembler(const ConstraintAssembler&) = delete;

    ConstraintAssembler& operator=(const ConstraintAssembler&) = delete;

    static const Variable<std::string>& GetImpositionVariable() noexcept
    {
        return IDENTIFIER;
    }

    typename TSparse::MatrixType mRelationMatrix;

    typename TSparse::VectorType mConstraintGapVector;

    std::string mName;
}; // class ConstraintImposition


} // namespace Kratos
