---
type: C++ Function
title: ComputeArea2DGeometry
resource: utilities/integration_utilities.h#L107-L119
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/Det2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/quadrilateral_2d_4/Quadrilateral2D4/Area
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_2d_8/Quadrilateral2D8/Area
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_2d_9/Quadrilateral2D9/Area
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_6/Triangle2D6/Area
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_integration_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double ComputeArea2DGeometry(const Geometry<TPointType>& rGeometry)`

# Calls

- [GetDefaultIntegrationMethod](../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [Det2](../../../functions/utilities/math_utils/Det2.md)
- [Weight](../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [Area](../../../functions/geometries/quadrilateral_2d_4/Quadrilateral2D4/Area.md)
- [Area](../../../functions/geometries/quadrilateral_2d_8/Quadrilateral2D8/Area.md)
- [Area](../../../functions/geometries/quadrilateral_2d_9/Quadrilateral2D9/Area.md)
- [Area](../../../functions/geometries/triangle_2d_6/Triangle2D6/Area.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_integration_utilities/KRATOS_TEST_CASE_IN_SUITE.md)