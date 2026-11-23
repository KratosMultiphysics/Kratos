---
type: C++ Method
title: SwapElementNode
resource: processes/mesh_node_collapsing_process.cpp#L133-L140
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CollapseNode
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MeshNodeCollapsingProcess::SwapElementNode(Element& rElement, Node const& rThisNode, Node::Pointer pCoarseNode)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CollapseNode](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CollapseNode.md)