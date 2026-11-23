---
type: C++ Method
title: GetDistance
resource: spatial_containers/spatial_search_result.h#L127-L129
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_nodal_distance_to_skin_process/CalculateNodalDistanceToSkinProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchNearestInRadius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchNearestInCell
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/spatial_containers/test_geometrical_objects_bins/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/spatial_containers/test_spatial_search_result/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double GetDistance() const`

# Called by

- [Execute](../../../../functions/processes/calculate_nodal_distance_to_skin_process/CalculateNodalDistanceToSkinProcess/Execute.md)
- [SearchNearestInRadius](../../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchNearestInRadius.md)
- [SearchNearestInCell](../../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchNearestInCell.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/spatial_containers/test_geometrical_objects_bins/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/spatial_containers/test_spatial_search_result/KRATOS_TEST_CASE_IN_SUITE.md)