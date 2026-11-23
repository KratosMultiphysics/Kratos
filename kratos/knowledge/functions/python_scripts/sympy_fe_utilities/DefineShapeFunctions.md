---
type: Python Function
title: DefineShapeFunctions
resource: python_scripts/sympy_fe_utilities.py#L66-L102
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/sympy_fe_utilities/DefineMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/DefineVector
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testDefineShapeFunctions1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testDefineTestFunctions2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testDefineTestFunctions3
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def DefineShapeFunctions(nnodes, dim, impose_partion_of_unity=False, shape_functions_name='N', first_derivatives_name='DN', second_derivatives_name=None):`

# Calls

- [DefineMatrix](../../../functions/python_scripts/sympy_fe_utilities/DefineMatrix.md)
- [DefineVector](../../../functions/python_scripts/sympy_fe_utilities/DefineVector.md)

# Called by

- [testDefineShapeFunctions1](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testDefineShapeFunctions1.md)
- [testDefineTestFunctions2](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testDefineTestFunctions2.md)
- [testDefineTestFunctions3](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testDefineTestFunctions3.md)