---
type: Python Function
title: DefineVector
resource: python_scripts/sympy_fe_utilities.py#L56-L64
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/sympy_fe_utilities/DefineShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testVectorOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testOutputScalar_CollectingFactors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testAuxiliaryOutputCollectionFactorsWithPower
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def DefineVector( name, m):`

# Called by

- [DefineShapeFunctions](../../../functions/python_scripts/sympy_fe_utilities/DefineShapeFunctions.md)
- [testVectorOutput](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testVectorOutput.md)
- [testOutputScalar_CollectingFactors](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testOutputScalar_CollectingFactors.md)
- [testAuxiliaryOutputCollectionFactorsWithPower](../../../functions/tests/test_sympy_fe_utilities/TestSympyFEUtilities/testAuxiliaryOutputCollectionFactorsWithPower.md)