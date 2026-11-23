---
type: Python Function
title: _AuxiliaryOutputCollectionFactors
resource: python_scripts/sympy_fe_utilities.py#L742-L769
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/sympy_fe_utilities/OutputSymbolicVariableDeclaration
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/sympy_fe_utilities/OutputMatrix_CollectingFactors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/OutputVector_CollectingFactors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/OutputScalar_CollectingFactors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testAuxiliaryOutputCollectionFactorsWithPower
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _AuxiliaryOutputCollectionFactors(A, name, language, indentation_level, optimizations, replace_indices, assignment_op, output_func):`

# Calls

- [OutputSymbolicVariableDeclaration](../../../functions/python_scripts/sympy_fe_utilities/OutputSymbolicVariableDeclaration.md)

# Called by

- [OutputMatrix_CollectingFactors](../../../functions/python_scripts/sympy_fe_utilities/OutputMatrix_CollectingFactors.md)
- [OutputVector_CollectingFactors](../../../functions/python_scripts/sympy_fe_utilities/OutputVector_CollectingFactors.md)
- [OutputScalar_CollectingFactors](../../../functions/python_scripts/sympy_fe_utilities/OutputScalar_CollectingFactors.md)
- [testAuxiliaryOutputCollectionFactorsWithPower](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testAuxiliaryOutputCollectionFactorsWithPower.md)