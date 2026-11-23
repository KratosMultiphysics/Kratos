---
type: C++ Method
title: Execute
resource: processes/mesh_node_collapsing_process.cpp#L36-L43
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/RemoveElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MeshNodeCollapsingProcess::Execute()`

# Calls

- [RemoveElements](../../../../functions/sources/model_part/ModelPart/RemoveElements.md)
- [RemoveNodes](../../../../functions/sources/model_part/ModelPart/RemoveNodes.md)