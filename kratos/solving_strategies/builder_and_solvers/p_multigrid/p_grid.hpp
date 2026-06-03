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
#include "solving_strategies/builder_and_solvers/p_multigrid/constraint_assembler.hpp" // ConstraintAssembler
#include "solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp" // DiagonalScaling
#include "solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities.hpp" // detail::DofData
#include "linear_solvers/linear_solver.h" // LinearSolver, Reorderer
#include "includes/dof.h" // Dof

// System includes
#include <memory> // std::shared_ptr
#include <optional> // std::optional
#include <vector> // std::vector


namespace Kratos {


/** @brief Coarse grid level of @ref PMultigridBuilderAndSolver.
 *  @details Settings:
 *           @code
 *           {
 *               "max_depth" : 0,
 *               "precision" : "double",
 *               "verbosity" : 1,
 *               "constraint_imposition_settings" : {
 *                   "method" : "augmented_lagrange",
 *                   "max_iterations" : 1
 *               }
 *           }
 *           @endcode
 *           - @p "max_depth" Maximum number of coarse grid levels. If set to 0 (default), there are no
 *                            coarse grids and no multigrid preconditioning is done. This also means that
 *                            that the top-level grid uses the linear solver instead of the smoother.
 *           - @p "precision" Floating point precision of system matrices and vectors. Either @p "double" (default)
 *                            or @p "single". The potential benefit of using single precision systems is reducing
 *                            memory consumtion on accelerator hardware (e.g.: VRAM on GPUs). This is viable because
 *                            coarse corrections in a p-multigrid preconditioner need not have high accuracy.
 *           - @p "verbosity" Level of verbosity. Every level includes lower verbosity levels and
 *                            adds new events to report:
 *              - @p 0 No messages, not even in case of failure.
 *              - @p 1 Warnings and failure messages.
 *              - @p 2 Aggregated status reports.
 *              - @p 3 Per-iteration status reports.
 *              - @p 4 Output system matrices and vectors.
 *          - @p "constraint_imposition_settings" : Settings of the top-level constraint assembler.
 *                                                  Refer to @ref ConstraintAssemblerFactory for more
 *                                                  information.
 */
template <class TSparse, class TDense>
class PGrid
{
public:
    using LinearSolverType = LinearSolver<TSparse,TDense,Reorderer<TSparse,TDense>>;

    using IndirectDofSet = PointerVectorSet<Dof<typename TDense::DataType>>;

    PGrid();

    PGrid(Parameters Settings,
          Parameters SmootherSettings,
          Parameters LeafSolverSettings,
          Parameters DiagonalScalingSettings);

    PGrid(PGrid&&) noexcept = default;

    PGrid& operator=(PGrid&&) noexcept = default;

    template <class TParentSparse>
    void MakeLhsTopology(ModelPart& rModelPart,
                         const typename TParentSparse::MatrixType& rParentLhs,
                         const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler,
                         const IndirectDofSet& rParentDofSet);

    template <bool AssembleLHS,
              bool AssembleRHS,
              class TParentSparse>
    void Assemble(ModelPart& rModelPart,
                  const typename TParentSparse::MatrixType* pParentLhs,
                  const typename TParentSparse::VectorType* pParentRhs,
                  const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler,
                  IndirectDofSet& rParentDofSet);

    void ApplyDirichletConditions(typename IndirectDofSet::const_iterator itParentDofBegin,
                                  typename IndirectDofSet::const_iterator itParentDofEnd);

    void ApplyConstraints();

    template <class TParentSparse>
    void Initialize(ModelPart& rModelPart,
                    const typename TParentSparse::MatrixType& rParentLhs,
                    const typename TParentSparse::VectorType& rParentSolution,
                    const typename TParentSparse::VectorType& rParentRhs);

    template <class TParentSparse>
    bool ApplyCoarseCorrection(typename TParentSparse::VectorType& rParentSolution,
                               const typename TParentSparse::VectorType& rParentResidual,
                               const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler,
                               PMGStatusStream& rStream);

    template <class TParentSparse>
    void Finalize(ModelPart& rModelPart,
                  const typename TParentSparse::MatrixType& rParentLhs,
                  const typename TParentSparse::VectorType& rParentSolution,
                  const typename TParentSparse::VectorType& rParentRhs);

    /// @brief Compute the residual in the coarse grid's independent space.
    /// @tparam TParentSparse Sparse space type of the parent grid.
    /// @param rCoarseIndependentResidual Output vector. Residual in the coarse grid's independent space.
    /// @param rFineIndependentRhs Residual vector in the fine grid's independent space.
    /// @param rParentConstraintAssembler Constraint assembler of the fine grid.
    template <class TParentSparse>
    void Restrict(typename TSparse::VectorType& rCoarseIndependentResidual,
                  const typename TParentSparse::VectorType& rFineIndependentResidual,
                  const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler) const;

    /// @brief Compute the solution in the fine grid's independent space.
    /// @tparam TParentSparse Space space type of the parent grid.
    /// @param rFineIndependentSolution Output vector. Solution vector in the fine grid's independent space.
    /// @param rCoarseIndependentSolution Solution vector in the coarse grid's independent space.
    /// @param rParentConstraintAssembler Constraint assembler of the fine grid.
    template <class TParentSparse>
    void Prolong(typename TParentSparse::VectorType& rFineIndependentSolution,
                 const typename TSparse::VectorType& rCoarseIndependentSolution,
                 const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler) const;

    void Clear();

    static Parameters GetDefaultParameters();

    std::optional<const PGrid*> GetChild() const
    {
        return mMaybeChild.has_value() ? mMaybeChild.value().get() : std::optional<const PGrid*>();
    }

    const typename TSparse::VectorType& GetSolution() const
    {
        return mSolution;
    }

private:
    PGrid(Parameters Settings,
          const unsigned CurrentDepth,
          Parameters SmootherSettings,
          Parameters LeafSolverSettings,
          Parameters DiagonalScalingSettings);

    PGrid(const PGrid&) = delete;

    PGrid& operator=(const PGrid&) = delete;

    void ExecuteMultigridLoop(PMGStatusStream& rStream,
                              PMGStatusStream::Report& rReport);

    void ExecuteConstraintLoop(PMGStatusStream& rStream,
                               PMGStatusStream::Report& rReport);

    typename TSparse::MatrixType mRestrictionOperator;

    typename TSparse::MatrixType mProlongationOperator;

    typename TSparse::MatrixType mLhs;

    typename TSparse::VectorType mSolution;

    typename TSparse::VectorType mRhs;

    VariablesList::Pointer mpVariableList;

    /// @details Array of @ref Dof "DoFs" unique to the current grid level.
    ///          DoFs need a pointer to a @ref NodalData object, which is why
    ///          pairs of @ref NodalData and @ref Dof are stored instead of just
    ///          DoFs.
    std::vector<detail::DofData> mDofSet;

    IndirectDofSet mIndirectDofSet;

    std::vector<std::size_t> mDofMap;

    std::shared_ptr<ConstraintAssembler<TSparse,TDense>> mpConstraintAssembler;

    std::unique_ptr<Scaling> mpDiagonalScaling;

    typename LinearSolverType::Pointer mpSolver;

    std::optional<std::unique_ptr<PGrid>> mMaybeChild;

    int mVerbosity;

    unsigned mDepth;
}; // class PGrid


} // namespace Kratos
