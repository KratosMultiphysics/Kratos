---
type: Python Function
title: OutputScalar
resource: python_scripts/sympy_fe_utilities.py#L613-L636
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
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testScalarOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def OutputScalar(scalar_expression, name, language, indentation_level=0, replace_indices=True, assignment_op="="):`

# Calls

- [_Indentation](../../../functions/python_scripts/sympy_fe_utilities/_Indentation.md)
- [_Suffix](../../../functions/python_scripts/sympy_fe_utilities/_Suffix.md)
- [_CodeGen](../../../functions/python_scripts/sympy_fe_utilities/_CodeGen.md)
- [_ReplaceIndices](../../../functions/python_scripts/sympy_fe_utilities/_ReplaceIndices.md)

# Called by

- [testScalarOutput](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testScalarOutput.md)