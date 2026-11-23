---
type: Python Method
title: test_matrix_interface
resource: tests/test_kratos_parameters.py#L889-L966
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_matrix_interface(self): # Read and check Matrices from a Parameters-Object`

# Calls

- [IsMatrix](../../../../functions/sources/kratos_parameters/Parameters/IsMatrix.md)
- [assertMatrixAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)
- [AddEmptyValue](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [SetMatrix](../../../../functions/sources/kratos_parameters/Parameters/SetMatrix.md)