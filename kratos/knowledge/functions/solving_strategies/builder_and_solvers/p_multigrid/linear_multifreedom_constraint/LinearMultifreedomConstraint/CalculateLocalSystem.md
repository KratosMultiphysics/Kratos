---
type: C++ Method
title: CalculateLocalSystem
resource: solving_strategies/builder_and_solvers/p_multigrid/linear_multifreedom_constraint.cpp#L49-L55
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void LinearMultifreedomConstraint::CalculateLocalSystem(MatrixType& rRelationMatrix, VectorType& rConstraintGaps, const ProcessInfo&) const`

# Calls

- [GetData](../../../../../../functions/sources/gid_io/GidIOBase/GetData.md)