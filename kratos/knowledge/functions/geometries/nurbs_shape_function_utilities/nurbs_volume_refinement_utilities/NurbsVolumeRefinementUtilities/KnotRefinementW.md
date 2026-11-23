---
type: C++ Method
title: KnotRefinementW
resource: geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities.h#L334-L462
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfKnotsW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetLowerSpan
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/geometries/test_nurbs_volume_knot_refinement/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void KnotRefinementW(NurbsVolumeGeometryType& rGeometry, std::vector<double>& rKnotsWToInsert, PointerVector<NodeType>& rPointsRefined, Vector& rKnotsWRefined )`

# Calls

- [NumberOfKnotsW](../../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfKnotsW.md)
- [NumberOfControlPointsW](../../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW.md)
- [GetLowerSpan](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetLowerSpan.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_volume_knot_refinement/KRATOS_TEST_CASE_IN_SUITE.md)