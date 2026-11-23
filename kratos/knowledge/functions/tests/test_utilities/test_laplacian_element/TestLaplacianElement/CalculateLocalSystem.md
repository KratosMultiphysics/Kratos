---
type: C++ Method
title: CalculateLocalSystem
resource: tests/test_utilities/test_laplacian_element.cpp#L167-L217
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/JacobianOnInitialConfiguration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/ShapeFunctionsGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TestLaplacianElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [JacobianOnInitialConfiguration](../../../../../functions/utilities/geometry_utilities/JacobianOnInitialConfiguration.md)
- [ShapeFunctionsGradients](../../../../../functions/utilities/geometry_utilities/ShapeFunctionsGradients.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [prod](../../../../../functions/includes/amatrix_interface/prod.md)