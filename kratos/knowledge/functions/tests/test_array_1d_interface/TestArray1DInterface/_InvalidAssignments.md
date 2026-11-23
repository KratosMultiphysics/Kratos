---
type: Python Method
title: _InvalidAssignments
resource: tests/test_array_1d_interface.py#L128-L145
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_array_1d_interface/TestArray1DInterface/test_SetNodalArrayValueFromPython_WrongInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_array_1d_interface/TestArray1DInterface/test_SetElementalArrayValueFromPython_WrongInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_array_1d_interface/TestArray1DInterface/test_SetConditionalArrayValueFromPython_WrongInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_array_1d_interface/TestArray1DInterface/test_SetPropertyArrayValueFromPython_WrongInput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _InvalidAssignments(self, tested_object): # Note: the expected behavior is to obtain a pybind-defined TypeError complaining about "incompatible function arguments". # However, in debug compilations ublas may detect size inconsistencies first and throw a ValueError (bad argument) instead. # Forbidden implicit conversion from wrong-size Array`

# Called by

- [test_SetNodalArrayValueFromPython_WrongInput](../../../../functions/tests/test_array_1d_interface/TestArray1DInterface/test_SetNodalArrayValueFromPython_WrongInput.md)
- [test_SetElementalArrayValueFromPython_WrongInput](../../../../functions/tests/test_array_1d_interface/TestArray1DInterface/test_SetElementalArrayValueFromPython_WrongInput.md)
- [test_SetConditionalArrayValueFromPython_WrongInput](../../../../functions/tests/test_array_1d_interface/TestArray1DInterface/test_SetConditionalArrayValueFromPython_WrongInput.md)
- [test_SetPropertyArrayValueFromPython_WrongInput](../../../../functions/tests/test_array_1d_interface/TestArray1DInterface/test_SetPropertyArrayValueFromPython_WrongInput.md)