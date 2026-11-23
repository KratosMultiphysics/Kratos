---
type: C++ Method
title: InitializeConstraintIteration
resource: solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler.cpp#L390-L413
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/GetPenaltyFactor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AugmentedLagrangeConstraintAssembler<TSparse,TDense>::InitializeConstraintIteration(typename TSparse::MatrixType& rLhs, typename TSparse::VectorType& rSolution, typename TSparse::VectorType& rRhs, [[maybe_unused]] typename Base::DofSet::iterator itDofBegin, [[maybe_unused]] typename Base::DofSet::iterator itDofEnd)`

# Calls

- [GetPenaltyFactor](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/GetPenaltyFactor.md)