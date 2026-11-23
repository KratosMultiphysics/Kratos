---
type: Python Method
title: test_mm_vector_io
resource: tests/test_matrix_market_interface.py#L22-L33
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/matrix_market_interface/ReadMatrixMarketVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_mm_vector_io(self):`

# Calls

- [ReadMatrixMarketVector](../../../../functions/sources/matrix_market_interface/ReadMatrixMarketVector.md)
- [assertVectorAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)
- [DeleteFileIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)