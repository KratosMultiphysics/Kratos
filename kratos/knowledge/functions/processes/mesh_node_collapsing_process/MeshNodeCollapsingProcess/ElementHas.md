---
type: C++ Method
title: ElementHas
resource: processes/mesh_node_collapsing_process.cpp#L125-L131
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
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CalculateMinQualityOfNeighbourElements
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool MeshNodeCollapsingProcess::ElementHas(Element& rElement, Node const& rCoarseNode)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CollapseNode](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CollapseNode.md)
- [CalculateMinQualityOfNeighbourElements](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CalculateMinQualityOfNeighbourElements.md)