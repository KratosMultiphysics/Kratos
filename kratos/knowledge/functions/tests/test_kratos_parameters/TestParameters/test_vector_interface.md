---
type: Python Method
title: test_vector_interface
resource: tests/test_kratos_parameters.py#L818-L887
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_vector_interface(self): # Read and check Vectors from a Parameters-Object`

# Calls

- [IsVector](../../../../functions/sources/kratos_parameters/Parameters/IsVector.md)
- [assertVectorAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)
- [AddEmptyValue](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [SetVector](../../../../functions/sources/kratos_parameters/Parameters/SetVector.md)