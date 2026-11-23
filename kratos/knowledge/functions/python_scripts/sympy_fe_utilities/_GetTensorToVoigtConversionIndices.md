---
type: Python Function
title: _GetTensorToVoigtConversionIndices
resource: python_scripts/sympy_fe_utilities.py#L592-L611
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/sympy_fe_utilities/ConvertTensorToVoigtMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/_AuxCheckToMatrixResults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _GetTensorToVoigtConversionIndices(dim):`

# Called by

- [ConvertTensorToVoigtMatrix](../../../functions/python_scripts/sympy_fe_utilities/ConvertTensorToVoigtMatrix.md)
- [_AuxCheckToMatrixResults](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/_AuxCheckToMatrixResults.md)