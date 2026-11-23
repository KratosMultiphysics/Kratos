---
type: C++ Method
title: FillNodalConnectivitiesFromGeometryBlockInList
resource: sources/model_part_io.cpp#L1121-L1185
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
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalGraphFromEntitiesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::FillNodalConnectivitiesFromGeometryBlockInList( ConnectivitiesContainerType& rNodalConnectivities, std::unordered_set<SizeType>& rGeometriesIds)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReorderedGeometryId](../../../../functions/sources/model_part_io/ModelPartIO/ReorderedGeometryId.md)

# Called by

- [ReadNodalGraphFromEntitiesList](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraphFromEntitiesList.md)