---
type: Python Function
title: OutputMatrix
resource: python_scripts/sympy_fe_utilities.py#L665-L693
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/sympy_fe_utilities/_Indentation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/_Suffix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/_CodeGen
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/_ReplaceIndices
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testMatrixOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testListOfMatricesOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def OutputMatrix(matrix_expression, name, language, indentation_level=0, replace_indices=True, assignment_op="="):`

# Calls

- [_Indentation](../../../functions/python_scripts/sympy_fe_utilities/_Indentation.md)
- [_Suffix](../../../functions/python_scripts/sympy_fe_utilities/_Suffix.md)
- [_CodeGen](../../../functions/python_scripts/sympy_fe_utilities/_CodeGen.md)
- [_ReplaceIndices](../../../functions/python_scripts/sympy_fe_utilities/_ReplaceIndices.md)

# Called by

- [testMatrixOutput](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testMatrixOutput.md)
- [testListOfMatricesOutput](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testListOfMatricesOutput.md)