---
type: Python Method
title: test_scalar_prod
resource: tests/test_matrix_interface.py#L100-L110
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

`def test_scalar_prod(self):`

# Calls

- [fill](../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [assertMatrixAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)