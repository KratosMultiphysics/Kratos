---
type: C++ Method
title: ReadGeometriesConnectivitiesBlock
resource: sources/model_part_io.cpp#L3354-L3408
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReorderedGeometryId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPartIO::SizeType ModelPartIO::ReadGeometriesConnectivitiesBlock(ConnectivitiesContainerType& rThisConnectivities)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReorderedGeometryId](../../../../functions/sources/model_part_io/ModelPartIO/ReorderedGeometryId.md)

# Called by

- [ReadGeometriesConnectivities](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivities.md)