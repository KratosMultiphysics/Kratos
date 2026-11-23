---
type: Python Function
title: to_csr
resource: python_scripts/scipy_conversion_tools.py#L8-L10
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_scipy_conversion_tools/TestScipyConversionTools/test_scipy_conversion_tools
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_matrix_assembly
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def to_csr(A):`

# Called by

- [test_scipy_conversion_tools](../../../functions/tests/test_scipy_conversion_tools/TestScipyConversionTools/test_scipy_conversion_tools.md)
- [test_matrix_assembly](../../../functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_matrix_assembly.md)