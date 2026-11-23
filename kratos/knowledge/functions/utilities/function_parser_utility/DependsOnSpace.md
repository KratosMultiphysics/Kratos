---
type: C++ Function
title: DependsOnSpace
resource: utilities/function_parser_utility.h#L116-L119
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<TEntity, THistorical>/InternalAssignValueVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<TEntity, THistorical>/InternalAssignValueScalar
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/set_initial_state_process/components_to_functions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_apply_function_to_nodes_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_function_parser_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_GenericFunctionUtility1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_GenericFunctionUtility2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_ApplyFunctionToNodesUtility
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool DependsOnSpace()`

# Called by

- [InternalAssignValueVector](../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<TEntity, THistorical>/InternalAssignValueVector.md)
- [InternalAssignValueScalar](../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<TEntity, THistorical>/InternalAssignValueScalar.md)
- [__init__](../../../functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/__init__.md)
- [ExecuteInitializeSolutionStep](../../../functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/ExecuteInitializeSolutionStep.md)
- [components_to_functions](../../../functions/python_scripts/set_initial_state_process/components_to_functions.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_apply_function_to_nodes_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_function_parser_utility/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_GenericFunctionUtility1](../../../functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_GenericFunctionUtility1.md)
- [test_GenericFunctionUtility2](../../../functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_GenericFunctionUtility2.md)
- [test_ApplyFunctionToNodesUtility](../../../functions/tests/test_function_parser_utility/TestGenericFunctionUtility/test_ApplyFunctionToNodesUtility.md)