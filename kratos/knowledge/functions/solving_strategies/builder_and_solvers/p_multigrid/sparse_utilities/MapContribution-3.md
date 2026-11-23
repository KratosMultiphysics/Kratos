---
type: C++ Function
title: MapContribution
resource: solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp#L360-L378
generated:
  by: okf-rs/0.7.0
---

# Signature

`void MapContribution(typename TSparse::MatrixType& rLhs, typename TSparse::VectorType& rRhs, const typename TDense::MatrixType& rLhsContribution, const typename TDense::VectorType& rRhsContribution, const Element::EquationIdVectorType& rEquationIds, LockObject* pLockBegin) noexcept`