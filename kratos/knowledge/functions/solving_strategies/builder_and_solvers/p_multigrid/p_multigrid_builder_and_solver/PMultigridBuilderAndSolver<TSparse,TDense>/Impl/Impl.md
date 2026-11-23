---
type: C++ Method
title: Impl
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L85-L96
generated:
  by: okf-rs/0.7.0
---

# Signature

`Impl(Interface* pInterface) : mpInterface(pInterface), mpMaybeModelPart(), mpConstraintAssembler(), mMaybeHierarchy(), mpRootGridSolver(), mpDiagonalScaling(), mMaxIterations(0), mTolerance(std::numeric_limits<typename TSparse::DataType>::max()), mMaxDepth(-1), mVerbosity(0)`