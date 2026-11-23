---
type: Python Method
title: _AuxCheckToVoigtResults
resource: tests/test_sympy_fe_utilities.py#L128-L134
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/sympy_fe_utilities/_GetVoigtToTensorConversionIndices
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testStrainToVoigt2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testStrainToVoigt3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testMatrixToVoigt2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testMatrixToVoigt3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _AuxCheckToVoigtResults(self, strain_size, factor_calculator, reference, voigt_output):`

# Calls

- [_GetVoigtToTensorConversionIndices](../../../../functions/python_scripts/sympy_fe_utilities/_GetVoigtToTensorConversionIndices.md)

# Called by

- [testStrainToVoigt2D](../../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testStrainToVoigt2D.md)
- [testStrainToVoigt3D](../../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testStrainToVoigt3D.md)
- [testMatrixToVoigt2D](../../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testMatrixToVoigt2D.md)
- [testMatrixToVoigt3D](../../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testMatrixToVoigt3D.md)