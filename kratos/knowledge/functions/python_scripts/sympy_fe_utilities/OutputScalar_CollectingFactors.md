---
type: Python Function
title: OutputScalar_CollectingFactors
resource: python_scripts/sympy_fe_utilities.py#L812-L825
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/sympy_fe_utilities/_AuxiliaryOutputCollectionFactors
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testOutputScalar_CollectingFactors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def OutputScalar_CollectingFactors(A, name, language, indentation_level=0, optimizations='basic', replace_indices=True, assignment_op="="):`

# Calls

- [_AuxiliaryOutputCollectionFactors](../../../functions/python_scripts/sympy_fe_utilities/_AuxiliaryOutputCollectionFactors.md)

# Called by

- [testOutputScalar_CollectingFactors](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testOutputScalar_CollectingFactors.md)