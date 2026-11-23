---
type: Python Method
title: test_mm_matrix_io
resource: tests/test_matrix_market_interface.py#L7-L20
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/matrix_market_interface/ReadMatrixMarketMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_mm_matrix_io(self):`

# Calls

- [ReadMatrixMarketMatrix](../../../../functions/sources/matrix_market_interface/ReadMatrixMarketMatrix.md)
- [assertMatrixAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)
- [DeleteFileIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)