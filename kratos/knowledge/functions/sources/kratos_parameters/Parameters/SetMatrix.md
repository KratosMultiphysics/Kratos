---
type: C++ Method
title: SetMatrix
resource: sources/kratos_parameters.cpp#L888-L904
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/csv_points_output_process/Factory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/line_output_process/LineOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_set_methods
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_matrix_interface
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::SetMatrix(const Matrix& rValue)`

# Called by

- [Factory](../../../../functions/python_scripts/csv_points_output_process/Factory.md)
- [__init__](../../../../functions/python_scripts/line_output_process/LineOutputProcess/__init__.md)
- [AddMatrix](../../../../functions/sources/kratos_parameters/Parameters/AddMatrix.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_set_methods](../../../../functions/tests/test_kratos_parameters/TestParameters/test_set_methods.md)
- [test_matrix_interface](../../../../functions/tests/test_kratos_parameters/TestParameters/test_matrix_interface.md)