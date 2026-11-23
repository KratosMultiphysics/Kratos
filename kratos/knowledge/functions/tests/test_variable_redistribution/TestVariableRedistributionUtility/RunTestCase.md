---
type: Python Method
title: RunTestCase
resource: tests/test_variable_redistribution.py#L159-L190
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_SetUpProblem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_CheckDoubleResults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_CheckArrayResults
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_linear_function
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_sharp_corners
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_quadratic
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_vector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def RunTestCase(self,inteface_check,set_reference,reference_variable,intermediate_variable,result_variable): # Read mdpa and set test model part`

# Calls

- [_SetUpProblem](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_SetUpProblem.md)
- [_PrintOutput](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_PrintOutput.md)
- [_CheckDoubleResults](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_CheckDoubleResults.md)
- [_CheckArrayResults](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_CheckArrayResults.md)

# Called by

- [test_linear_function](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_linear_function.md)
- [test_sharp_corners](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_sharp_corners.md)
- [test_quadratic](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_quadratic.md)
- [test_vector](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_vector.md)