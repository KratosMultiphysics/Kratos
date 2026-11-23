---
type: C++ Method
title: IsElementInsideBounds
resource: modeler/operation/find_contacts_in_skin_model_part.cpp#L244-L247
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetGridSize
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetCellNeihgbourColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetCellNeihgbourColorsConnectedByFace
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool FindContactsInSkinModelPart::IsElementInsideBounds(const CartesianMeshColors& rColors, int i, int j, int k) const`

# Calls

- [GetGridSize](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetGridSize.md)

# Called by

- [GetCellNeihgbourColors](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetCellNeihgbourColors.md)
- [GetCellNeihgbourColorsConnectedByFace](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetCellNeihgbourColorsConnectedByFace.md)