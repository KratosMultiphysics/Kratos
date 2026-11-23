---
type: C++ Method
title: DivideNodalDataBlock
resource: sources/model_part_io.cpp#L5171-L5218
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/WriteInAllFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideFlagVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideDofVariableData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::DivideNodalDataBlock(OutputFilesContainerType& OutputFiles, PartitionIndicesContainerType const& NodesAllPartitions)`

# Calls

- [WriteInAllFiles](../../../../functions/sources/model_part_io/ModelPartIO/WriteInAllFiles.md)
- [DivideFlagVariableData](../../../../functions/sources/model_part_io/ModelPartIO/DivideFlagVariableData.md)
- [DivideDofVariableData](../../../../functions/sources/model_part_io/ModelPartIO/DivideDofVariableData.md)

# Called by

- [DivideInputToPartitionsImpl](../../../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl.md)