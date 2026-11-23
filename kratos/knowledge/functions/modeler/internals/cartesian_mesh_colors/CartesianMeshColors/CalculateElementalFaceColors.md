---
type: C++ Method
title: CalculateElementalFaceColors
resource: modeler/internals/cartesian_mesh_colors.h#L388-L553
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/CollapseIntersectionPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_ray/CartesianRay/MarkIntersectedIntervals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/coloring/color_cell_faces/ColorCellFaces/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateElementalFaceColors(array_1d< std::size_t, 3 > const& MinRayPosition, array_1d< std::size_t, 3 > const& MaxRayPosition, int FaceColor, int OutsideColor, int VolumeColor)`

# Calls

- [CollapseIntersectionPoints](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/CollapseIntersectionPoints.md)
- [MarkIntersectedIntervals](../../../../../functions/modeler/internals/cartesian_ray/CartesianRay/MarkIntersectedIntervals.md)
- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)

# Called by

- [Apply](../../../../../functions/modeler/coloring/color_cell_faces/ColorCellFaces/Apply.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE.md)