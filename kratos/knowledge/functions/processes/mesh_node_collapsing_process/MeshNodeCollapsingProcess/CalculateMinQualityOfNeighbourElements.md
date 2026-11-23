---
type: C++ Method
title: CalculateMinQualityOfNeighbourElements
resource: processes/mesh_node_collapsing_process.cpp#L113-L123
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/ElementHas
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CalculateQualityIfNodeCollapses
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double MeshNodeCollapsingProcess::CalculateMinQualityOfNeighbourElements(Node& rThisNode, Node const& rCoarseNode)`

# Calls

- [ElementHas](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/ElementHas.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CalculateQualityIfNodeCollapses](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CalculateQualityIfNodeCollapses.md)