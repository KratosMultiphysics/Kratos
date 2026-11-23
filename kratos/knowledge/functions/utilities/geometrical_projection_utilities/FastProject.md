---
type: C++ Function
title: FastProject
resource: utilities/geometrical_projection_utilities.h#L131-L147
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/triangle_3d_3/IsInside
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_triangle_3d_3/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometrical_projection_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, false>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, false>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, false, 4>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, false, 3>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, true>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, true>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, true, 4>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, true, 3>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/ProjectedIsInside
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline TPointClass3 FastProject( const TPointClass1& rPointOrigin, const TPointClass2& rPointToProject, const array_1d<double,3>& rNormal, double& rDistance )`

# Called by

- [IsInside](../../../functions/geometries/triangle_3d_3/IsInside.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_triangle_3d_3/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_geometrical_projection_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, false>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, false>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, false, 4>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, false, 3>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, true>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, true>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, true, 4>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, true, 3>/GetExactIntegration.md)
- [ProjectedIsInside](../../../functions/utilities/geometry_utilities/GeometryUtils/ProjectedIsInside.md)