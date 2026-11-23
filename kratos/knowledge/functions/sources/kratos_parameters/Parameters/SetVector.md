---
type: C++ Method
title: SetVector
resource: sources/kratos_parameters.cpp#L873-L883
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/tests/test_linear_system/TestLinearSystem/test_LinearSystemWithLinearOperator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/multiple_points_output_process/MultiplePointsOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_initial_time/TestInitialTime/_SetInitialTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_set_methods
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_vector_interface
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::SetVector(const Vector& rValue)`

# Called by

- [test_LinearSystemWithLinearOperator](../../../../functions/future/tests/test_linear_system/TestLinearSystem/test_LinearSystemWithLinearOperator.md)
- [__init__](../../../../functions/python_scripts/multiple_points_output_process/MultiplePointsOutputProcess/__init__.md)
- [AddVector](../../../../functions/sources/kratos_parameters/Parameters/AddVector.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [_SetInitialTime](../../../../functions/tests/test_initial_time/TestInitialTime/_SetInitialTime.md)
- [test_set_methods](../../../../functions/tests/test_kratos_parameters/TestParameters/test_set_methods.md)
- [test_vector_interface](../../../../functions/tests/test_kratos_parameters/TestParameters/test_vector_interface.md)