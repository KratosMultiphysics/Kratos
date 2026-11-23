---
type: Python Method
title: _CheckDoubleResults
resource: tests/test_variable_redistribution.py#L237-L241
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_quad_elements_quadratic_scalar_field
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/RunTestCase
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_nodal_area
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _CheckDoubleResults(self, model_part, reference_variable,result_variable):`

# Called by

- [test_quad_elements_quadratic_scalar_field](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_quad_elements_quadratic_scalar_field.md)
- [RunTestCase](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/RunTestCase.md)
- [test_nodal_area](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_nodal_area.md)