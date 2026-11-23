---
type: Python Method
title: __CheckValues
resource: tests/test_tensor_adaptors.py#L1057-L1065
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestVariableTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CheckValues(self, value_1, value_2):`

# Calls

- [assertVectorAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)
- [assertMatrixAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)

# Called by

- [__TestVariableTensorAdaptor](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/__TestVariableTensorAdaptor.md)