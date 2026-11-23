---
type: Python Method
title: _AuxCheckToMatrixResults
resource: tests/test_sympy_fe_utilities.py#L136-L141
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/sympy_fe_utilities/_GetTensorToVoigtConversionIndices
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testVoigtToMatrix2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testVoigtToMatrix3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testStrainToMatrix2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testStrainToMatrix3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _AuxCheckToMatrixResults(self, dimension, factor_calculator, reference, matrix_output):`

# Calls

- [_GetTensorToVoigtConversionIndices](../../../../functions/python_scripts/sympy_fe_utilities/_GetTensorToVoigtConversionIndices.md)

# Called by

- [testVoigtToMatrix2D](../../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testVoigtToMatrix2D.md)
- [testVoigtToMatrix3D](../../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testVoigtToMatrix3D.md)
- [testStrainToMatrix2D](../../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testStrainToMatrix2D.md)
- [testStrainToMatrix3D](../../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testStrainToMatrix3D.md)