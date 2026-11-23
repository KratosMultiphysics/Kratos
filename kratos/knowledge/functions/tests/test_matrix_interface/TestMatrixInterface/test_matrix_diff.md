---
type: Python Method
title: test_matrix_diff
resource: tests/test_matrix_interface.py#L86-L98
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_matrix_diff(self):`

# Calls

- [fill](../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [assertMatrixAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)