---
type: C++ Method
title: IsStringArray
resource: sources/kratos_parameters.cpp#L584-L597
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesToCheckpointList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_is_methods
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Parameters::IsStringArray() const`

# Called by

- [ReplaceDeprecatedVariableName](../../../../functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName.md)
- [__GetStagesToCheckpointList](../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesToCheckpointList.md)
- [GetStringArray](../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_is_methods](../../../../functions/tests/test_kratos_parameters/TestParameters/test_is_methods.md)