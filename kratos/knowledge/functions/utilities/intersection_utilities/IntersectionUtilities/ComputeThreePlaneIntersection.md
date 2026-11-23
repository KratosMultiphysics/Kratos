---
type: C++ Method
title: ComputeThreePlaneIntersection
resource: utilities/intersection_utilities.cpp#L72-L150
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/intersection_utilities/IntersectionUtilities/ComputeTwoPlaneIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_intersection_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int IntersectionUtilities::ComputeThreePlaneIntersection( const array_1d<double, 3>& rPlanePoint1, const array_1d<double, 3>& rNormal1, const array_1d<double, 3>& rPlanePoint2, const array_1d<double, 3>& rNormal2, const array_1d<double, 3>& rPlanePoint3, const array_1d<double, 3>& rNormal3, array_1d<double, 3>& rPoint1, array_1d<double, 3>& rVector1, array_1d<double, 3>& rPoint2, array_1d<double, 3>& rVector2, const double Tolerance )`

# Calls

- [ComputeTwoPlaneIntersection](../../../../functions/utilities/intersection_utilities/IntersectionUtilities/ComputeTwoPlaneIntersection.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_intersection_utilities/KRATOS_TEST_CASE_IN_SUITE.md)