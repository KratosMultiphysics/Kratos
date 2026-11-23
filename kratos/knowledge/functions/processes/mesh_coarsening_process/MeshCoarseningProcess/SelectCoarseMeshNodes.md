---
type: C++ Method
title: SelectCoarseMeshNodes
resource: processes/mesh_coarsening_process.cpp#L50-L82
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/mesh_coarsening_process/MeshCoarseningProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MeshCoarseningProcess::SelectCoarseMeshNodes()`

# Calls

- [IsNot](../../../../functions/containers/flags/IsNot.md)

# Called by

- [Execute](../../../../functions/processes/mesh_coarsening_process/MeshCoarseningProcess/Execute.md)