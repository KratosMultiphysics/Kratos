---
type: Python Function
title: _ReplaceIndices
resource: python_scripts/sympy_fe_utilities.py#L541-L565
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/sympy_fe_utilities/OutputScalar
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/OutputVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/OutputMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/OutputSymbolicVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/OutputSymbolicVariableDeclaration
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _ReplaceIndices(language, expression):`

# Called by

- [OutputScalar](../../../functions/python_scripts/sympy_fe_utilities/OutputScalar.md)
- [OutputVector](../../../functions/python_scripts/sympy_fe_utilities/OutputVector.md)
- [OutputMatrix](../../../functions/python_scripts/sympy_fe_utilities/OutputMatrix.md)
- [OutputSymbolicVariable](../../../functions/python_scripts/sympy_fe_utilities/OutputSymbolicVariable.md)
- [OutputSymbolicVariableDeclaration](../../../functions/python_scripts/sympy_fe_utilities/OutputSymbolicVariableDeclaration.md)