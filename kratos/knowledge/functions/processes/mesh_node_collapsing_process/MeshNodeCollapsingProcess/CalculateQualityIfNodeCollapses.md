---
type: C++ Method
title: CalculateQualityIfNodeCollapses
resource: processes/mesh_node_collapsing_process.cpp#L105-L111
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CalculateMinQualityOfNeighbourElements
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CollapseNode
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double MeshNodeCollapsingProcess::CalculateQualityIfNodeCollapses(Node& rThisNode, Node const& rCoarseNode)`

# Calls

- [CalculateMinQualityOfNeighbourElements](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CalculateMinQualityOfNeighbourElements.md)

# Called by

- [CollapseNode](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CollapseNode.md)