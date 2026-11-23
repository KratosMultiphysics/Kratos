---
type: C++ Method
title: GetElementalColor
resource: modeler/internals/cartesian_mesh_colors.h#L129-L131
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/coloring/color_and_intersect_with_cells_in_touch/ColorAndIntersectWithCellsInTouch/ApplyColorToCellsInTouchWithGeometryAndComputeIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_cells_in_touch/ColorCellsInTouch/ApplyColorToCellsInTouchWithGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_cells_with_center_in_sphere_around_nodes/ColorCellsWithCenterInSphereArounNodes/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ApplyColorToConnectedCellsInTouchWithGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ColorConnectedCellsToThisCell
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_outer_faces_of_cells_with_colors/ColorOuterFacesOfCellsWithColors/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_outer_faces_of_cells_with_colors/ColorOuterFacesOfCellsWithColors/ApplyColorIfOuterFace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/generate_elements_with_cell_color/GenerateElementsWithCellColor/Generate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/generate_quadratic_elements_with_cell_color/GenerateQuadraticElementsWithCellColor/Generate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/Generate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/generate_tetrahedral_elements_with_cell_color/GenerateTetrahedralElementsWithCellColor/Generate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalRayColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalFaceColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalFaceColorsBetweenColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/PrintAllVoxelColorsToAsciiFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/FindConditionContact
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetCellNeihgbourColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetCellNeihgbourColorsConnectedByFace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double GetElementalColor(std::size_t I, std::size_t J, std::size_t K) const`

# Called by

- [ApplyColorToCellsInTouchWithGeometryAndComputeIntersection](../../../../../functions/modeler/coloring/color_and_intersect_with_cells_in_touch/ColorAndIntersectWithCellsInTouch/ApplyColorToCellsInTouchWithGeometryAndComputeIntersection.md)
- [ApplyColorToCellsInTouchWithGeometry](../../../../../functions/modeler/coloring/color_cells_in_touch/ColorCellsInTouch/ApplyColorToCellsInTouchWithGeometry.md)
- [Apply](../../../../../functions/modeler/coloring/color_cells_with_center_in_sphere_around_nodes/ColorCellsWithCenterInSphereArounNodes/Apply.md)
- [ApplyColorToConnectedCellsInTouchWithGeometry](../../../../../functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ApplyColorToConnectedCellsInTouchWithGeometry.md)
- [ColorConnectedCellsToThisCell](../../../../../functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ColorConnectedCellsToThisCell.md)
- [Apply](../../../../../functions/modeler/coloring/color_outer_faces_of_cells_with_colors/ColorOuterFacesOfCellsWithColors/Apply.md)
- [ApplyColorIfOuterFace](../../../../../functions/modeler/coloring/color_outer_faces_of_cells_with_colors/ColorOuterFacesOfCellsWithColors/ApplyColorIfOuterFace.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_elements_with_cell_color/GenerateElementsWithCellColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_quadratic_elements_with_cell_color/GenerateQuadraticElementsWithCellColor/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/Generate.md)
- [Generate](../../../../../functions/modeler/entity_generation/generate_tetrahedral_elements_with_cell_color/GenerateTetrahedralElementsWithCellColor/Generate.md)
- [CalculateElementalRayColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalRayColors.md)
- [CalculateElementalFaceColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalFaceColors.md)
- [CalculateElementalFaceColorsBetweenColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateElementalFaceColorsBetweenColors.md)
- [PrintAllVoxelColorsToAsciiFile](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/PrintAllVoxelColorsToAsciiFile.md)
- [FindConditionContact](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/FindConditionContact.md)
- [GetCellNeihgbourColors](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetCellNeihgbourColors.md)
- [GetCellNeihgbourColorsConnectedByFace](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetCellNeihgbourColorsConnectedByFace.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE.md)