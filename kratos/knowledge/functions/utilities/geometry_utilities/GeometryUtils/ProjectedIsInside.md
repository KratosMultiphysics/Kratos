---
type: C++ Method
title: ProjectedIsInside
resource: utilities/geometry_utilities.cpp#L694-L723
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometrical_projection_utilities/FastProjectOnLine2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/FastProject
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryUtils::ProjectedIsInside( const GeometryType& rGeometry, const GeometryType::CoordinatesArrayType& rPointGlobalCoordinates, GeometryType::CoordinatesArrayType& rResult, const double Tolerance )`

# Calls

- [FastProjectOnLine2D](../../../../functions/utilities/geometrical_projection_utilities/FastProjectOnLine2D.md)
- [FastProject](../../../../functions/utilities/geometrical_projection_utilities/FastProject.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE.md)