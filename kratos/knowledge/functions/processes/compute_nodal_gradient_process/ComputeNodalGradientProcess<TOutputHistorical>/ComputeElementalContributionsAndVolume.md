---
type: C++ Method
title: ComputeElementalContributionsAndVolume
resource: processes/compute_nodal_gradient_process.cpp#L119-L217
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
  - target: functions/includes/model_part/FullName
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
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeNodalGradientProcess<TOutputHistorical>::ComputeElementalContributionsAndVolume()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)
- [GetDefaultIntegrationMethod](../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [JacobianOnInitialConfiguration](../../../../functions/utilities/geometry_utilities/JacobianOnInitialConfiguration.md)
- [GeneralizedInvertMatrix](../../../../functions/utilities/math_utils/GeneralizedInvertMatrix.md)
- [ShapeFunctionsGradients](../../../../functions/utilities/geometry_utilities/ShapeFunctionsGradients.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)
- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [Execute](../../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/Execute.md)