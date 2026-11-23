---
type: C++ Method
title: ComputeLineBoxIntersection
resource: utilities/intersection_utilities.cpp#L362-L391
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/intersection_utilities/GetLineBoxIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/intersection_utilities/InBox
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/line_3d_2/Line3D2/HasIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_intersection_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int IntersectionUtilities::ComputeLineBoxIntersection( const array_1d<double,3>& rBoxPoint0, const array_1d<double,3>& rBoxPoint1, const array_1d<double,3>& rLinePoint0, const array_1d<double,3>& rLinePoint1)`

# Calls

- [GetLineBoxIntersection](../../../../functions/utilities/intersection_utilities/GetLineBoxIntersection.md)
- [InBox](../../../../functions/utilities/intersection_utilities/InBox.md)

# Called by

- [HasIntersection](../../../../functions/geometries/line_3d_2/Line3D2/HasIntersection.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_intersection_utilities/KRATOS_TEST_CASE_IN_SUITE.md)