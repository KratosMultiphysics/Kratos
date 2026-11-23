---
type: C++ Method
title: GetCenterOfElement
resource: modeler/internals/cartesian_mesh_colors.h#L161-L163
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/coloring/color_cells_with_center_in_sphere_around_nodes/ColorCellsWithCenterInSphereArounNodes/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetClosestContactColorFromNeighbours
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Point GetCenterOfElement(std::size_t I, std::size_t J, std::size_t K) const`

# Called by

- [Apply](../../../../../functions/modeler/coloring/color_cells_with_center_in_sphere_around_nodes/ColorCellsWithCenterInSphereArounNodes/Apply.md)
- [GetClosestContactColorFromNeighbours](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetClosestContactColorFromNeighbours.md)