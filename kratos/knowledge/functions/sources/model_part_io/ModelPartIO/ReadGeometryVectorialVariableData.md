---
type: C++ Method
title: ReadGeometryVectorialVariableData
resource: sources/model_part_io.cpp#L3226-L3255
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
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometryDataBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadGeometryVectorialVariableData(GeometryContainerType& rThisGeometries, const TVariableType& rVariable, TDataType Dummy)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReorderedGeometryId](../../../../functions/sources/model_part_io/ModelPartIO/ReorderedGeometryId.md)

# Called by

- [ReadGeometryDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometryDataBlock.md)