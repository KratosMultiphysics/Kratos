---
type: C++ Method
title: IsBool
resource: sources/kratos_parameters.cpp#L560-L563
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/assign_vector_variable_process/AssignVectorVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesToCheckpointList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_kratos_change_parameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_is_methods
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Parameters::IsBool() const`

# Called by

- [__init__](../../../../functions/python_scripts/assign_vector_variable_process/AssignVectorVariableProcess/__init__.md)
- [ReplaceDeprecatedVariableName](../../../../functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName.md)
- [__GetStagesToCheckpointList](../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesToCheckpointList.md)
- [__PrepareOutputSettings](../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_kratos_change_parameters](../../../../functions/tests/test_kratos_parameters/TestParameters/test_kratos_change_parameters.md)
- [test_is_methods](../../../../functions/tests/test_kratos_parameters/TestParameters/test_is_methods.md)