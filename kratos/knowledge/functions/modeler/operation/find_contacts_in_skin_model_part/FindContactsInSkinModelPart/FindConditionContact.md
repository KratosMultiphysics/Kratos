---
type: C++ Method
title: FindConditionContact
resource: modeler/operation/find_contacts_in_skin_model_part.cpp#L81-L159
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CalculateOuterMinMaxNodePositions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetCellNeihgbourColorsConnectedByFace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetClosestContactColorFromNeighbours
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetCellNeihgbourColors
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindContactsInSkinModelPart::FindConditionContact( Condition& rCondition, const int CellColor, ContactContainer& rContactContainer) const`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CalculateOuterMinMaxNodePositions](../../../../../functions/modeler/internals/cartesian_mesh_colors/CalculateOuterMinMaxNodePositions.md)
- [GetElementalColor](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/GetElementalColor.md)
- [GetCellNeihgbourColorsConnectedByFace](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetCellNeihgbourColorsConnectedByFace.md)
- [GetClosestContactColorFromNeighbours](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetClosestContactColorFromNeighbours.md)
- [GetCellNeihgbourColors](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetCellNeihgbourColors.md)

# Called by

- [Execute](../../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/Execute.md)