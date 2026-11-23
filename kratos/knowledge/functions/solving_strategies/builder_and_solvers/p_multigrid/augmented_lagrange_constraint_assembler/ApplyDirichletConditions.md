---
type: C++ Function
title: ApplyDirichletConditions
resource: solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler.cpp#L32-L125
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/dof/Dof/IsFree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FixDof
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyDirichletConditions(typename TSparse::MatrixType& rRelationMatrix, typename TSparse::VectorType& rConstraintGaps, TItDof itDofBegin, [[maybe_unused]] TItDof itDofEnd, const std::string& rConstraintAssemblerName, int Verbosity)`

# Calls

- [IsFree](../../../../../functions/includes/dof/Dof/IsFree.md)
- [FixDof](../../../../../functions/includes/dof/Dof/FixDof.md)