---
type: C++ Method
title: LinearMultifreedomConstraint
resource: solving_strategies/builder_and_solvers/p_multigrid/linear_multifreedom_constraint.cpp#L21-L46
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`LinearMultifreedomConstraint::LinearMultifreedomConstraint(const IndexType Id, DofPointerVectorType&& rDofs, const std::vector<std::size_t>& rConstraintLabels, const MatrixType& rRelationMatrix, const VectorType& rConstraintGaps) : MultifreedomConstraint(Id, std::move(rDofs), rConstraintLabels)`