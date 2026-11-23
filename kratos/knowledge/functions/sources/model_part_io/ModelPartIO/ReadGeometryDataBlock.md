---
type: C++ Method
title: ReadGeometryDataBlock
resource: sources/model_part_io.cpp#L3160-L3190
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometryScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometryVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadInitialValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadGeometryDataBlock(GeometryContainerType& rThisGeometries)`

# Calls

- [ReadGeometryScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometryScalarVariableData.md)
- [ReadGeometryVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometryVectorialVariableData.md)

# Called by

- [ReadInitialValues](../../../../functions/sources/model_part_io/ModelPartIO/ReadInitialValues.md)