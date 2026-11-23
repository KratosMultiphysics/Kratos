---
type: Python Method
title: ExecuteInitializeSolutionStep
resource: python_scripts/assign_scalar_variable_process.py#L104-L129
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/interval_utility/IntervalUtility/IsInInterval
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/DependsOnSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/apply_function_to_nodes_utility/ApplyFunctionToNodesUtility/ApplyFunction
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ExecuteInitializeSolutionStep(self):`

# Calls

- [IsInInterval](../../../../functions/utilities/interval_utility/IntervalUtility/IsInInterval.md)
- [DependsOnSpace](../../../../functions/utilities/function_parser_utility/DependsOnSpace.md)
- [ApplyFunction](../../../../functions/utilities/apply_function_to_nodes_utility/ApplyFunctionToNodesUtility/ApplyFunction.md)