---
type: Python Function
title: comparator
resource: tests/test_serializer.py#L172-L199
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
  - target: functions/tests/test_serializer/TestSerializer/test_DataOnly
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def comparator(model_1: KratosMultiphysics.Model, model_2: KratosMultiphysics.Model, factor: int):`

# Calls

- [assertVectorAlmostEqual](../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)
- [assertMatrixAlmostEqual](../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)

# Called by

- [test_DataOnly](../../../functions/tests/test_serializer/TestSerializer/test_DataOnly.md)