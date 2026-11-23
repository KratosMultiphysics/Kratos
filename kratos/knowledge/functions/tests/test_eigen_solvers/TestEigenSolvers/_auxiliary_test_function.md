---
type: Python Method
title: _auxiliary_test_function
resource: tests/test_eigen_solvers.py#L19-L48
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/matrix_market_interface/ReadMatrixMarketMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _auxiliary_test_function(self, settings, matrix_name="auxiliar_files_for_python_unittest/sparse_matrix_files/A.mm", eigen_value_estimated = "lowest"):`

# Calls

- [ReadMatrixMarketMatrix](../../../../functions/sources/matrix_market_interface/ReadMatrixMarketMatrix.md)