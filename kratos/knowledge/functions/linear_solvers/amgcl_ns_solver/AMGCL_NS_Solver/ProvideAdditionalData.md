---
type: C++ Method
title: ProvideAdditionalData
resource: linear_solvers/amgcl_ns_solver.h#L317-L419
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ProvideAdditionalData ( SparseMatrixType& rA, VectorType& rX, VectorType& rB, typename ModelPart::DofsArrayType& rDofSet, ModelPart& rModelPart ) override`

# Calls

- [GetVariable](../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)