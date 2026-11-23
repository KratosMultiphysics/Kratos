---
type: C++ Function
title: FastProjectOnGeometry
resource: utilities/geometrical_projection_utilities.h#L159-L176
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometrical_projection_utilities/FastProjectDirection
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_geometrical_projection_utilities/TestFastProjectOnGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double FastProjectOnGeometry(const TGeometryType& rGeom, const Point& rPointToProject, Point& rPointProjected, const SizeType EchoLevel = 0)`

# Calls

- [FastProjectDirection](../../../functions/utilities/geometrical_projection_utilities/FastProjectDirection.md)

# Called by

- [TestFastProjectOnGeometry](../../../functions/tests/cpp_tests/utilities/test_geometrical_projection_utilities/TestFastProjectOnGeometry.md)