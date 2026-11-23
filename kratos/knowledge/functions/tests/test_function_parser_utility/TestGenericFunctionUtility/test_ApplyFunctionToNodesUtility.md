---
type: Python Method
title: test_ApplyFunctionToNodesUtility
resource: tests/test_function_parser_utility.py#L96-L127
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/function_parser_utility/DependsOnSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/BasicGenericFunctionUtility/FunctionBody
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/apply_function_to_nodes_utility/ApplyFunctionToNodesUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/apply_function_to_nodes_utility/ApplyFunctionToNodesUtility/ApplyFunction
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_ApplyFunctionToNodesUtility(self):`

# Calls

- [DependsOnSpace](../../../../functions/utilities/function_parser_utility/DependsOnSpace.md)
- [FunctionBody](../../../../functions/utilities/function_parser_utility/BasicGenericFunctionUtility/FunctionBody.md)
- [ApplyFunctionToNodesUtility](../../../../functions/utilities/apply_function_to_nodes_utility/ApplyFunctionToNodesUtility.md)
- [ApplyFunction](../../../../functions/utilities/apply_function_to_nodes_utility/ApplyFunctionToNodesUtility/ApplyFunction.md)