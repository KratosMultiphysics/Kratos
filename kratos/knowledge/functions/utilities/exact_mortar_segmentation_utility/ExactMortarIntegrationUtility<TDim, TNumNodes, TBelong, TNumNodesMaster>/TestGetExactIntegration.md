---
type: C++ Method
title: TestGetExactIntegration
resource: utilities/exact_mortar_segmentation_utility.cpp#L1073-L1098
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_line_exact_integration_1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_line_exact_integration_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_line_exact_integration_3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/_test_triangle_exact_integration_1
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::TestGetExactIntegration( Condition::Pointer pSlaveCond, Condition::Pointer pMasterCond, Matrix& rCustomSolution )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [test_line_exact_integration_1](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_line_exact_integration_1.md)
- [test_line_exact_integration_2](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_line_exact_integration_2.md)
- [test_line_exact_integration_3](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_line_exact_integration_3.md)
- [_test_triangle_exact_integration_1](../../../../functions/tests/test_exact_integration/TestExactIntegration/_test_triangle_exact_integration_1.md)