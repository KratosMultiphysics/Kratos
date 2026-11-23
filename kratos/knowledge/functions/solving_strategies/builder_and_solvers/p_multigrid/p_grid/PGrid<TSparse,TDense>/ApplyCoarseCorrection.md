---
type: C++ Method
title: ApplyCoarseCorrection
resource: solving_strategies/builder_and_solvers/p_multigrid/p_grid.cpp#L460-L494
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`bool PGrid<TSparse,TDense>::ApplyCoarseCorrection(typename TParentSparse::VectorType& rParentSolution, const typename TParentSparse::VectorType& rParentRhs, const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler, PMGStatusStream& rStream)`