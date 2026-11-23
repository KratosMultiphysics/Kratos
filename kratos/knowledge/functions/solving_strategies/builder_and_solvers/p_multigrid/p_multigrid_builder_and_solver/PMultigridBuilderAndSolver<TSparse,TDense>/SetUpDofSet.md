---
type: C++ Method
title: SetUpDofSet
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L576-L672
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/HasReaction
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PMultigridBuilderAndSolver<TSparse,TDense>::SetUpDofSet(typename Interface::TSchemeType::Pointer pScheme, ModelPart& rModelPart)`

# Calls

- [fill](../../../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [GetVariable](../../../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [HasReaction](../../../../../../functions/includes/dof/Dof/HasReaction.md)