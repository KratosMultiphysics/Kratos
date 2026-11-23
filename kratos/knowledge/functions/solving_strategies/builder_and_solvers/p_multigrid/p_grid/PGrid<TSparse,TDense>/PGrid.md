---
type: C++ Method
title: PGrid
resource: solving_strategies/builder_and_solvers/p_multigrid/p_grid.cpp#L34-L121
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`PGrid<TSparse,TDense>::PGrid(Parameters Settings, const unsigned CurrentDepth, Parameters SmootherSettings, Parameters LeafSolverSettings, Parameters DiagonalScalingSettings) : mRestrictionOperator(), mProlongationOperator(), mLhs(), mSolution(), mRhs(), mDofSet(), mIndirectDofSet(), mDofMap(), mpConstraintAssembler(), mpDiagonalScaling(), mpSolver(), mMaybeChild(), mVerbosity(), mDepth(CurrentDepth)`