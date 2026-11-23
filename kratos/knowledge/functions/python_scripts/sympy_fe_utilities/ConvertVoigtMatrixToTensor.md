---
type: Python Function
title: ConvertVoigtMatrixToTensor
resource: python_scripts/sympy_fe_utilities.py#L236-L262
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/sympy_fe_utilities/_GetVoigtToTensorConversionIndices
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testConvertVoigtMatrixToTensor2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testConvertVoigtMatrixToTensor3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ConvertVoigtMatrixToTensor(C):`

# Calls

- [_GetVoigtToTensorConversionIndices](../../../functions/python_scripts/sympy_fe_utilities/_GetVoigtToTensorConversionIndices.md)

# Called by

- [testConvertVoigtMatrixToTensor2D](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testConvertVoigtMatrixToTensor2D.md)
- [testConvertVoigtMatrixToTensor3D](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testConvertVoigtMatrixToTensor3D.md)