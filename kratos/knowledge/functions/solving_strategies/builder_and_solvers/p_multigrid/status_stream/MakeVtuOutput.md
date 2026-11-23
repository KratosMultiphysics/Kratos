---
type: C++ Function
title: MakeVtuOutput
resource: solving_strategies/builder_and_solvers/p_multigrid/status_stream.cpp#L42-L65
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/AddVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Impl/WriteIntermediateState
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::unique_ptr<VtuOutput> MakeVtuOutput(ModelPart& rModelPart, const PointerVectorSet<Dof<double>>& rDofSet)`

# Calls

- [GetVariable](../../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [AddVariable](../../../../../functions/input_output/vtu_output/VtuOutput/AddVariable.md)

# Called by

- [WriteIntermediateState](../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Impl/WriteIntermediateState.md)