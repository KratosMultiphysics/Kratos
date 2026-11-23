---
type: C++ Function
title: MapEntityContribution
resource: solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp#L390-L445
generated:
  by: okf-rs/0.7.0
---

# Signature

`void MapEntityContribution(TEntity& rEntity, Scheme<TSparse,TDense>& rScheme, const ProcessInfo& rProcessInfo, Element::EquationIdVectorType& rEquationIndices, typename TDense::MatrixType* pLhsContribution, typename TDense::VectorType* pRhsContribution, typename TSparse::MatrixType* pLhs, typename TSparse::VectorType* pRhs, LockObject* pLockBegin)`