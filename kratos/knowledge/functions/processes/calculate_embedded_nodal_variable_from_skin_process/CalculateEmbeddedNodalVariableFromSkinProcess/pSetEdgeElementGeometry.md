---
type: C++ Method
title: pSetEdgeElementGeometry
resource: processes/calculate_embedded_nodal_variable_from_skin_process.h#L763-L772
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddIntersectedElementsModelPartElements
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Element::GeometryType::Pointer pSetEdgeElementGeometry( std::unordered_map<unsigned int, Node::Pointer>& rNodesMap, const Element::GeometryType &rCurrentEdgeGeometry, const std::pair<std::size_t, std::size_t> NewEdgeIds) const`

# Called by

- [AddIntersectedElementsModelPartElements](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddIntersectedElementsModelPartElements.md)