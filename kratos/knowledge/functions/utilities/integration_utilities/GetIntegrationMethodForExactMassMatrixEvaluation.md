---
type: C++ Function
title: GetIntegrationMethodForExactMassMatrixEvaluation
resource: utilities/integration_utilities.h#L42-L54
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/line_2d_3/Line2D3/Length
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_4/Line2D4/Length
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_5/Line2D5/Length
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_3d_3/Line3D3/Length
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static GeometryData::IntegrationMethod GetIntegrationMethodForExactMassMatrixEvaluation(const Geometry<TPointType>& rGeometry)`

# Calls

- [GetDefaultIntegrationMethod](../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)

# Called by

- [Length](../../../functions/geometries/line_2d_3/Line2D3/Length.md)
- [Length](../../../functions/geometries/line_2d_4/Line2D4/Length.md)
- [Length](../../../functions/geometries/line_2d_5/Line2D5/Length.md)
- [Length](../../../functions/geometries/line_3d_3/Line3D3/Length.md)