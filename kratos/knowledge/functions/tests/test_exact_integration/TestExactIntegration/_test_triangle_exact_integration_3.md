---
type: Python Method
title: _test_triangle_exact_integration_3
resource: tests/test_exact_integration.py#L266-L321
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_triangle_exact_integration_3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_triangle_exact_integration_3_delaunator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _test_triangle_exact_integration_3(self, use_delaunator = False):`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [test_triangle_exact_integration_3](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_triangle_exact_integration_3.md)
- [test_triangle_exact_integration_3_delaunator](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_triangle_exact_integration_3_delaunator.md)