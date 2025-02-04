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
#include "linear_solvers/linear_solver.h" // LinearSolver, Reorderer

// System includes
#include <memory> // std::shared_ptr
#include <optional> // std::optional


namespace Kratos {


template <class TSparse, class TDense>
class PGrid
{
public:
    using LinearSolverType = LinearSolver<TSparse,TDense,Reorderer<TSparse,TDense>>;

    using DofSet = std::vector<Dof<typename TDense::DataType>>;

    using IndirectDofSet = PointerVectorSet<typename DofSet::value_type>;

    PGrid();

    PGrid(Parameters Settings,
          Parameters SmootherSettings,
          Parameters LeafSolverSettings);

    PGrid(PGrid&&) noexcept = default;

    PGrid& operator=(PGrid&&) noexcept = default;

    template <class TParentSparse>
    void MakeLhsTopology(ModelPart& rModelPart,
                         const typename TParentSparse::MatrixType& rParentLhs,
                         const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler);

    template <bool AssembleLHS,
              bool AssembleRHS,
              class TParentSparse>
    void Assemble(const ModelPart& rModelPart,
                  const typename TParentSparse::MatrixType* pParentLhs,
                  const typename TParentSparse::VectorType* pParentRhs,
                  const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler);

    template <class TParentSparse>
    void Initialize(ModelPart& rModelPart,
                    const typename TParentSparse::MatrixType& rParentLhs,
                    const typename TParentSparse::VectorType& rParentSolution,
                    const typename TParentSparse::VectorType& rParentRhs);

    template <class TParentSparse>
    bool ApplyCoarseCorrection(typename TParentSparse::VectorType& rParentSolution,
                               const typename TParentSparse::VectorType& rParentResidual) const;

    template <class TParentSparse>
    void Finalize(ModelPart& rModelPart,
                  const typename TParentSparse::MatrixType& rParentLhs,
                  const typename TParentSparse::VectorType& rParentSolution,
                  const typename TParentSparse::VectorType& rParentRhs);

    void Clear();

    static Parameters GetDefaultParameters();

private:
    PGrid(Parameters Settings,
          const unsigned CurrentDepth,
          Parameters SmootherSettings,
          Parameters LeafSolverSettings);

    PGrid(const PGrid&) = delete;

    PGrid& operator=(const PGrid&) = delete;

    typename TSparse::MatrixType mRestrictionOperator;

    typename TSparse::MatrixType mProlongationOperator;

    typename TSparse::MatrixType mLhs;

    typename TSparse::DataType mSolution;

    typename TSparse::VectorType mRhs;

    DofSet mDofSet;

    IndirectDofSet mIndirectDofSet;

    std::shared_ptr<ConstraintAssembler<TSparse,TDense>> mpConstraintAssembler;

    typename LinearSolverType::Pointer mpSolver;

    std::optional<std::unique_ptr<PGrid>> mMaybeChild;

    int mVerbosity;

    unsigned mDepth;
}; // class PGrid


} // namespace Kratos
