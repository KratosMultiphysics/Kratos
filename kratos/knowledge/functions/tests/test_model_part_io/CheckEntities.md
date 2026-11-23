---
type: Python Function
title: CheckEntities
resource: tests/test_model_part_io.py#L394-L400
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_model_part_io/GetScalar
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_io/GetArray3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CheckEntities(entites, ids_to_check):`

# Calls

- [GetScalar](../../../functions/tests/test_model_part_io/GetScalar.md)
- [assertVectorAlmostEqual](../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)
- [GetArray3](../../../functions/tests/test_model_part_io/GetArray3.md)
- [assertMatrixAlmostEqual](../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)