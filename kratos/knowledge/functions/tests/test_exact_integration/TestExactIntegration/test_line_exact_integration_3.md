---
type: Python Method
title: test_line_exact_integration_3
resource: tests/test_exact_integration.py#L105-L152
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
---

# Signature

`def test_line_exact_integration_3(self):`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [TestGetExactIntegration](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestGetExactIntegration.md)