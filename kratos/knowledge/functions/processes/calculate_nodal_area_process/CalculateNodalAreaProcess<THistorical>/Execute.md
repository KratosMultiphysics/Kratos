---
type: C++ Method
title: Execute
resource: processes/calculate_nodal_area_process.cpp#L70-L139
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/SolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
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
  - target: functions/utilities/math_utils/GeneralizedDet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateNodalAreaProcess<THistorical>::Execute()`

# Calls

- [SolutionStepsDataHas](../../../../functions/includes/node/SolutionStepsDataHas.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [GetDefaultIntegrationMethod](../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [JacobianOnInitialConfiguration](../../../../functions/utilities/geometry_utilities/JacobianOnInitialConfiguration.md)
- [GeneralizedDet](../../../../functions/utilities/math_utils/GeneralizedDet.md)
- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)