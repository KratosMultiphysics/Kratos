---
type: C++ Method
title: FillNodalConnectivitiesFromGeometryBlock
resource: sources/model_part_io.cpp#L3644-L3703
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalGraph
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::FillNodalConnectivitiesFromGeometryBlock(ConnectivitiesContainerType& rNodalConnectivities)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadNodalGraph](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraph.md)