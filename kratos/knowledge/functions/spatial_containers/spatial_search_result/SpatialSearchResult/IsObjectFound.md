---
type: C++ Method
title: IsObjectFound
resource: spatial_containers/spatial_search_result.h#L162-L164
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchNearestInRadius
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

`bool IsObjectFound() const`

# Called by

- [Execute](../../../../functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/Execute.md)
- [SearchNearestInRadius](../../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchNearestInRadius.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/spatial_containers/test_geometrical_objects_bins/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/spatial_containers/test_spatial_search_result/KRATOS_TEST_CASE_IN_SUITE.md)