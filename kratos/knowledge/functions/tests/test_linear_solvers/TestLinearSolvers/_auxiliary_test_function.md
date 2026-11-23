---
type: Python Method
title: _auxiliary_test_function
resource: tests/test_linear_solvers.py#L18-L70
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/matrix_market_interface/ReadMatrixMarketMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_strategies_to_python/SetToZeroVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _auxiliary_test_function(self, settings, matrix_name="auxiliar_files_for_python_unittest/sparse_matrix_files/A.mm"):`

# Calls

- [ReadMatrixMarketMatrix](../../../../functions/sources/matrix_market_interface/ReadMatrixMarketMatrix.md)
- [SetToZeroVector](../../../../functions/python/add_strategies_to_python/SetToZeroVector.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)