---
type: C++ Method
title: CollapseNode
resource: processes/mesh_node_collapsing_process.cpp#L65-L103
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CalculateQualityIfNodeCollapses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/ElementHas
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/SwapElementNode
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CollapseNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MeshNodeCollapsingProcess::CollapseNode(Node& rThisNode)`

# Calls

- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CalculateQualityIfNodeCollapses](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CalculateQualityIfNodeCollapses.md)
- [ElementHas](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/ElementHas.md)
- [SwapElementNode](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/SwapElementNode.md)

# Called by

- [CollapseNodes](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CollapseNodes.md)