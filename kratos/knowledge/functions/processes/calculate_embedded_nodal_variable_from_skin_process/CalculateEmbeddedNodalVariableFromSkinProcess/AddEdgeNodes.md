---
type: C++ Method
title: AddEdgeNodes
resource: processes/calculate_embedded_nodal_variable_from_skin_process.h#L737-L761
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/pGetNodalSolutionStepVariablesList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/SetSolutionStepVariablesList
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddIntersectedElementsModelPartElements
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddEdgeNodes( const Geometry<Node> &rEdgeGeometry, ModelPart &rModelPart, std::unordered_map<unsigned int, Node::Pointer>& rNodesMap ) const`

# Calls

- [pGetNodalSolutionStepVariablesList](../../../../functions/includes/model_part/pGetNodalSolutionStepVariablesList.md)
- [SetSolutionStepVariablesList](../../../../functions/includes/node/SetSolutionStepVariablesList.md)

# Called by

- [AddIntersectedElementsModelPartElements](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddIntersectedElementsModelPartElements.md)