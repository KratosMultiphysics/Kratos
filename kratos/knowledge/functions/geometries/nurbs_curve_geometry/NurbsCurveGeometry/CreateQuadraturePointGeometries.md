---
type: C++ Method
title: CreateQuadraturePointGeometries
resource: geometries/nurbs_curve_geometry.h#L521-L591
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/GetFirstNonzeroControlPoint
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateQuadraturePointGeometries( GeometriesArrayType& rResultGeometries, IndexType NumberOfShapeFunctionDerivatives, const IntegrationPointsArrayType& rIntegrationPoints, IntegrationInfo& rIntegrationInfo) override`

# Calls

- [GetDefaultIntegrationMethod](../../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [GetFirstNonzeroControlPoint](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/GetFirstNonzeroControlPoint.md)