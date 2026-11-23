---
type: C++ Method
title: ReorderedGeometryId
resource: sources/model_part_io.cpp#L6718-L6723
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromGeometryBlockInList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometriesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometryScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometryVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivitiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideGeometriesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartGeometriesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartEntitiesIds
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPartIO::SizeType ModelPartIO::ReorderedGeometryId(ModelPartIO::SizeType GeometryId)`

# Called by

- [FillNodalConnectivitiesFromGeometryBlockInList](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromGeometryBlockInList.md)
- [ReadGeometriesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesBlock.md)
- [ReadGeometryScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometryScalarVariableData.md)
- [ReadGeometryVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometryVectorialVariableData.md)
- [ReadGeometriesConnectivitiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivitiesBlock.md)
- [DivideGeometriesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideGeometriesBlock.md)
- [DivideVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/DivideVectorialVariableData.md)
- [DivideScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/DivideScalarVariableData.md)
- [DivideSubModelPartGeometriesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartGeometriesBlock.md)
- [ReadSubModelPartEntitiesIds](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartEntitiesIds.md)