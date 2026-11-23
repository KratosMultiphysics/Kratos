---
type: Python Method
title: test_matrix_sum
resource: tests/test_matrix_interface.py#L72-L84
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

`def test_matrix_sum(self):`

# Calls

- [fill](../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [assertMatrixAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)