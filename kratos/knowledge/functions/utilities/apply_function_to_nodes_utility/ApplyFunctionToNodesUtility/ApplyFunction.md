---
type: C++ Method
title: ApplyFunction
resource: utilities/apply_function_to_nodes_utility.cpp#L25-L50
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/X0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Y0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Z0
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/ExecuteInitializeSolutionStep
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

`void ApplyFunctionToNodesUtility::ApplyFunction( const Variable<double>& rVariable, const double t )`

# Calls

- [X0](../../../../functions/includes/node/X0.md)
- [Y0](../../../../functions/includes/node/Y0.md)
- [Z0](../../../../functions/includes/node/Z0.md)

# Called by

- [ExecuteInitializeSolutionStep](../../../../functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/ExecuteInitializeSolutionStep.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_apply_function_to_nodes_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_ApplyFunctionToNodesUtility](../../../../functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_ApplyFunctionToNodesUtility.md)
- [test_ApplyFunctionToNodesUtilityTimeEvolutionTernary](../../../../functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_ApplyFunctionToNodesUtilityTimeEvolutionTernary.md)