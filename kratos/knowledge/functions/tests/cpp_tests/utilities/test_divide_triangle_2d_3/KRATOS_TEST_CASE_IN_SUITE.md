---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/utilities/test_divide_triangle_2d_3.cpp#L23-L171
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveSubdivisions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeSubdivisions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveInterfaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeInterfaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveInterfacesParentIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeInterfacesParentIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(DivideGeometryTriangle2D3Horizontal, KratosCoreFastSuite)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetPositiveSubdivisions](../../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveSubdivisions.md)
- [GetNegativeSubdivisions](../../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeSubdivisions.md)
- [GetPositiveInterfaces](../../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveInterfaces.md)
- [GetNegativeInterfaces](../../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeInterfaces.md)
- [GetPositiveInterfacesParentIds](../../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveInterfacesParentIds.md)
- [GetNegativeInterfacesParentIds](../../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetNegativeInterfacesParentIds.md)
- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)