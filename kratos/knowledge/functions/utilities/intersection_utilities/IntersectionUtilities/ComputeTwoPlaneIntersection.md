---
type: C++ Method
title: ComputeTwoPlaneIntersection
resource: utilities/intersection_utilities.cpp#L25-L67
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_intersection_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/intersection_utilities/IntersectionUtilities/ComputeThreePlaneIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int IntersectionUtilities::ComputeTwoPlaneIntersection( const array_1d<double, 3>& rPlanePoint1, const array_1d<double, 3>& rNormal1, const array_1d<double, 3>& rPlanePoint2, const array_1d<double, 3>& rNormal2, array_1d<double, 3>& rPoint, array_1d<double, 3>& rTangent, const double Tolerance )`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_intersection_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [ComputeThreePlaneIntersection](../../../../functions/utilities/intersection_utilities/IntersectionUtilities/ComputeThreePlaneIntersection.md)