---
type: C++ Method
title: Execute
resource: processes/compute_nodal_normal_divergence_process.cpp#L22-L121
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/JacobianOnInitialConfiguration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/GeneralizedInvertMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/ShapeFunctionsGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeNodalNormalDivergenceProcess<THistorical>::Execute()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [GetDefaultIntegrationMethod](../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [JacobianOnInitialConfiguration](../../../../functions/utilities/geometry_utilities/JacobianOnInitialConfiguration.md)
- [GeneralizedInvertMatrix](../../../../functions/utilities/math_utils/GeneralizedInvertMatrix.md)
- [ShapeFunctionsGradients](../../../../functions/utilities/geometry_utilities/ShapeFunctionsGradients.md)
- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)