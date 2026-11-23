---
type: C++ Method
title: DivideElementalDataBlock
resource: sources/model_part_io.cpp#L5398-L5443
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/WriteInAllFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::DivideElementalDataBlock(OutputFilesContainerType& OutputFiles, PartitionIndicesContainerType const& ElementsAllPartitions)`

# Calls

- [WriteInAllFiles](../../../../functions/sources/model_part_io/ModelPartIO/WriteInAllFiles.md)
- [DivideScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/DivideScalarVariableData.md)

# Called by

- [DivideInputToPartitionsImpl](../../../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl.md)