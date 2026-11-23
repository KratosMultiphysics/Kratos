---
type: C++ Method
title: GetCellNeihgbourColors
resource: modeler/operation/find_contacts_in_skin_model_part.cpp#L199-L218
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/IsElementInsideBounds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/FindConditionContact
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindContactsInSkinModelPart::GetCellNeihgbourColors( const array_1d<std::size_t,3>& rCellIndexes, std::vector<int>& rNeighbourColors, const CartesianMeshColors& rColors) const`

# Calls

- [IsElementInsideBounds](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/IsElementInsideBounds.md)
- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)

# Called by

- [FindConditionContact](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/FindConditionContact.md)