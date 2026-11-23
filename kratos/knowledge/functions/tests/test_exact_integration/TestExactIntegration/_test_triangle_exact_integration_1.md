---
type: Python Method
title: _test_triangle_exact_integration_1
resource: tests/test_exact_integration.py#L156-L212
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestGetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_triangle_exact_integration_1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_triangle_exact_integration_1_delanautor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _test_triangle_exact_integration_1(self, use_delaunator = False):`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [TestGetExactIntegration](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestGetExactIntegration.md)

# Called by

- [test_triangle_exact_integration_1](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_triangle_exact_integration_1.md)
- [test_triangle_exact_integration_1_delanautor](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_triangle_exact_integration_1_delanautor.md)