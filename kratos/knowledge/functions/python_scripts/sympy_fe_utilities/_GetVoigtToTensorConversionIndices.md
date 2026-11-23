---
type: Python Function
title: _GetVoigtToTensorConversionIndices
resource: python_scripts/sympy_fe_utilities.py#L567-L590
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/sympy_fe_utilities/ConvertVoigtMatrixToTensor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/SubstituteMatrixValueByVoigtSymbols
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/_AuxCheckToVoigtResults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _GetVoigtToTensorConversionIndices(strain_size):`

# Called by

- [ConvertVoigtMatrixToTensor](../../../functions/python_scripts/sympy_fe_utilities/ConvertVoigtMatrixToTensor.md)
- [SubstituteMatrixValueByVoigtSymbols](../../../functions/python_scripts/sympy_fe_utilities/SubstituteMatrixValueByVoigtSymbols.md)
- [_AuxCheckToVoigtResults](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/_AuxCheckToVoigtResults.md)