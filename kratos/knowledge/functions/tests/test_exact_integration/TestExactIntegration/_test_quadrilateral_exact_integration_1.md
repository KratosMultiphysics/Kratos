---
type: Python Method
title: _test_quadrilateral_exact_integration_1
resource: tests/test_exact_integration.py#L330-L372
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_quadrilateral_exact_integration_1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_quadrilateral_exact_integration_1_delaunator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _test_quadrilateral_exact_integration_1(self, use_delaunator = False):`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [test_quadrilateral_exact_integration_1](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_quadrilateral_exact_integration_1.md)
- [test_quadrilateral_exact_integration_1_delaunator](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_quadrilateral_exact_integration_1_delaunator.md)