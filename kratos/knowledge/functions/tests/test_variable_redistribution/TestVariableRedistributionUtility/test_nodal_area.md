---
type: Python Method
title: test_nodal_area
resource: tests/test_variable_redistribution.py#L192-L217
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_SetUpProblem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_CheckDoubleResults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_nodal_area(self):`

# Calls

- [_SetUpProblem](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_SetUpProblem.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [_PrintOutput](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_PrintOutput.md)
- [_CheckDoubleResults](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/_CheckDoubleResults.md)