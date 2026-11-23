---
type: C++ Method
title: Build
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h#L181-L314
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/GetDofPosition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Build( typename TSchemeType::Pointer pScheme, ModelPart& r_model_part, TSystemMatrixType& A, TSystemVectorType& b) override`

# Calls

- [BuiltinTimer](../../../../../functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer.md)
- [GetDofPosition](../../../../../functions/includes/node/GetDofPosition.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)