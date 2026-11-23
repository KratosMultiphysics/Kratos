---
type: C++ Function
title: CalculateOuterMinMaxNodePositions
resource: modeler/internals/cartesian_mesh_colors.h#L614-L628
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateNodePosition
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/coloring/color_and_intersect_with_cells_in_touch/ColorAndIntersectWithCellsInTouch/ApplyColorToCellsInTouchWithGeometryAndComputeIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_cells_in_touch/ColorCellsInTouch/ApplyColorToCellsInTouchWithGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ApplyColorToConnectedCellsInTouchWithGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/GenerateIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/FindConditionContact
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateOuterMinMaxNodePositions(TPointsContainerType const& Points, array_1d< std::size_t, 3 >& MinNodePosition, array_1d< std::size_t, 3 >& MaxNodePosition) const`

# Calls

- [CalculateMinMaxPoints](../../../../functions/modeler/internals/cartesian_mesh_colors/CalculateMinMaxPoints.md)
- [CalculateNodePosition](../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/CalculateNodePosition.md)

# Called by

- [ApplyColorToCellsInTouchWithGeometryAndComputeIntersection](../../../../functions/modeler/coloring/color_and_intersect_with_cells_in_touch/ColorAndIntersectWithCellsInTouch/ApplyColorToCellsInTouchWithGeometryAndComputeIntersection.md)
- [ApplyColorToCellsInTouchWithGeometry](../../../../functions/modeler/coloring/color_cells_in_touch/ColorCellsInTouch/ApplyColorToCellsInTouchWithGeometry.md)
- [ApplyColorToConnectedCellsInTouchWithGeometry](../../../../functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ApplyColorToConnectedCellsInTouchWithGeometry.md)
- [GenerateIntersection](../../../../functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/GenerateIntersection.md)
- [FindConditionContact](../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/FindConditionContact.md)