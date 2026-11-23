---
type: Python Function
title: OutputSymbolicVariableDeclaration
resource: python_scripts/sympy_fe_utilities.py#L713-L740
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/sympy_fe_utilities/_Indentation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/_CodeGen
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/_VariableDeclaration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/_Suffix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/_ReplaceIndices
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/sympy_fe_utilities/_AuxiliaryOutputCollectionFactors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testOutputSymbolicVariableDeclaration
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def OutputSymbolicVariableDeclaration(expression, name, language, indentation_level=0, replace_indices=True):`

# Calls

- [_Indentation](../../../functions/python_scripts/sympy_fe_utilities/_Indentation.md)
- [_CodeGen](../../../functions/python_scripts/sympy_fe_utilities/_CodeGen.md)
- [_VariableDeclaration](../../../functions/python_scripts/sympy_fe_utilities/_VariableDeclaration.md)
- [_Suffix](../../../functions/python_scripts/sympy_fe_utilities/_Suffix.md)
- [_ReplaceIndices](../../../functions/python_scripts/sympy_fe_utilities/_ReplaceIndices.md)

# Called by

- [_AuxiliaryOutputCollectionFactors](../../../functions/python_scripts/sympy_fe_utilities/_AuxiliaryOutputCollectionFactors.md)
- [testOutputSymbolicVariableDeclaration](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testOutputSymbolicVariableDeclaration.md)