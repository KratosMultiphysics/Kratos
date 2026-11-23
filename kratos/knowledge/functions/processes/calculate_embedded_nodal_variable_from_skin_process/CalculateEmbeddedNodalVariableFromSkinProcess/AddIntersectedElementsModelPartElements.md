---
type: C++ Method
title: AddIntersectedElementsModelPartElements
resource: processes/calculate_embedded_nodal_variable_from_skin_process.h#L494-L608
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/SetEdgePair
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddEdgeNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/pSetEdgeElementGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/GenerateIntersectedEdgesElementsModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddIntersectedElementsModelPartElements(ModelPart &rModelPart) const`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SetEdgePair](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/SetEdgePair.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [AddEdgeNodes](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddEdgeNodes.md)
- [pSetEdgeElementGeometry](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/pSetEdgeElementGeometry.md)

# Called by

- [GenerateIntersectedEdgesElementsModelPart](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/GenerateIntersectedEdgesElementsModelPart.md)