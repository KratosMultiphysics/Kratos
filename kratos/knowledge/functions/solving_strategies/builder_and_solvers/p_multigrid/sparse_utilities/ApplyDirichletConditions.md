---
type: C++ Function
title: ApplyDirichletConditions
resource: solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp#L452-L527
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ApplyDirichletConditions(typename TSparse::MatrixType& rLhs, typename TSparse::VectorType& rRhs, const TItRowDof itRowDofBegin, const TItRowDof itRowDofEnd, const TItColumnDof itColumnDofBegin, [[maybe_unused]] const TItColumnDof itColumnDofEnd, const typename TSparse::DataType DiagonalScaleFactor)`

# Calls

- [GetVariable](../../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)