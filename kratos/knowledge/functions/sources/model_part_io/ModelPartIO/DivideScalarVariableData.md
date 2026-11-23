---
type: C++ Method
title: DivideScalarVariableData
resource: sources/model_part_io.cpp#L5445-L5506
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReorderedGeometryId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/DivideElementalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideConditionalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideConstraintDataBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::DivideScalarVariableData(OutputFilesContainerType& OutputFiles, PartitionIndicesContainerType const& EntitiesPartitions, std::string BlockName)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReorderedConstraintId](../../../../functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId.md)
- [ReorderedGeometryId](../../../../functions/sources/model_part_io/ModelPartIO/ReorderedGeometryId.md)

# Called by

- [DivideElementalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideElementalDataBlock.md)
- [DivideConditionalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideConditionalDataBlock.md)
- [DivideConstraintDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideConstraintDataBlock.md)