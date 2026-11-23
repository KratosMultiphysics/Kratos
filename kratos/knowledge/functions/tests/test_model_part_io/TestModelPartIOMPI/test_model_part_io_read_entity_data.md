---
type: Python Method
title: test_model_part_io_read_entity_data
resource: tests/test_model_part_io.py#L350-L406
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

`def test_model_part_io_read_entity_data(self): # testing if the assignment of entity data works correctly in serial and MPI`

# Calls

- [GetScalar](../../../../functions/tests/test_model_part_io/GetScalar.md)
- [assertVectorAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)
- [GetArray3](../../../../functions/tests/test_model_part_io/GetArray3.md)
- [assertMatrixAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)