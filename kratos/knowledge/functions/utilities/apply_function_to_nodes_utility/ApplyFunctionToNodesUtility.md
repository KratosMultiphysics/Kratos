---
type: C++ Function
title: ApplyFunctionToNodesUtility
resource: utilities/apply_function_to_nodes_utility.h#L75-L81
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_apply_function_to_nodes_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_ApplyFunctionToNodesUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_ApplyFunctionToNodesUtilityTimeEvolutionTernary
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ApplyFunctionToNodesUtility( ModelPart::NodesContainerType& rNodes, GenericFunctionUtility::Pointer pFunction ): mrNodes(rNodes), mpFunction(pFunction)`

# Called by

- [__init__](../../../functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/__init__.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_apply_function_to_nodes_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_ApplyFunctionToNodesUtility](../../../functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_ApplyFunctionToNodesUtility.md)
- [test_ApplyFunctionToNodesUtilityTimeEvolutionTernary](../../../functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_ApplyFunctionToNodesUtilityTimeEvolutionTernary.md)