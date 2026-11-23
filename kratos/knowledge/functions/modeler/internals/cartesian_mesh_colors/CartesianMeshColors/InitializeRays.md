---
type: C++ Method
title: InitializeRays
resource: modeler/internals/cartesian_mesh_colors.h#L174-L211
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/coloring/color_cell_faces/ColorCellFaces/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_cells_with_inside_center/ColorCellsWithInsideCenter/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InitializeRays(array_1d< std::size_t, 3 > const& MinRayPosition, array_1d< std::size_t, 3 > const& MaxRayPosition, const std::string& EntititesToColor)`

# Called by

- [Apply](../../../../../functions/modeler/coloring/color_cell_faces/ColorCellFaces/Apply.md)
- [Apply](../../../../../functions/modeler/coloring/color_cells_with_inside_center/ColorCellsWithInsideCenter/Apply.md)
- [Execute](../../../../../functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/Execute.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE.md)