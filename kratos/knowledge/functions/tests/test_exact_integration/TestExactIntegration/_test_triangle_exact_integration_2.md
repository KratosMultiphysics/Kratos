---
type: Python Method
title: _test_triangle_exact_integration_2
resource: tests/test_exact_integration.py#L220-L258
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_triangle_exact_integration_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_triangle_exact_integration_2_delanautor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _test_triangle_exact_integration_2(self, use_delaunator = False):`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [test_triangle_exact_integration_2](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_triangle_exact_integration_2.md)
- [test_triangle_exact_integration_2_delanautor](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_triangle_exact_integration_2_delanautor.md)