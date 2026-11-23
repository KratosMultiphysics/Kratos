---
type: C++ Method
title: MultifreedomConstraint
resource: solving_strategies/builder_and_solvers/p_multigrid/multifreedom_constraint.cpp#L24-L34
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`MultifreedomConstraint::MultifreedomConstraint(const IndexType Id, DofArray&& rDofs, const std::vector<std::size_t>& rConstraintLabels) : MasterSlaveConstraint(Id), mDofs(std::move(rDofs))`