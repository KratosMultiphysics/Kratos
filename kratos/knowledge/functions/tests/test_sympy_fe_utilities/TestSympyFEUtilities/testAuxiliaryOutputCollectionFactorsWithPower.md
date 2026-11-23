---
type: Python Method
title: testAuxiliaryOutputCollectionFactorsWithPower
resource: tests/test_sympy_fe_utilities.py#L292-L307
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/sympy_fe_utilities/DefineVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sympy_fe_utilities/_AuxiliaryOutputCollectionFactors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def testAuxiliaryOutputCollectionFactorsWithPower(self): # Note: create_expand_pow_optimization is not expanding expressions # when there is a sum as shown in the last term of the test`

# Calls

- [DefineVector](../../../../functions/python_scripts/sympy_fe_utilities/DefineVector.md)
- [_AuxiliaryOutputCollectionFactors](../../../../functions/python_scripts/sympy_fe_utilities/_AuxiliaryOutputCollectionFactors.md)