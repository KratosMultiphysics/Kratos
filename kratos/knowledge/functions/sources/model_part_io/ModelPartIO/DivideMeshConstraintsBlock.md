---
type: C++ Method
title: DivideMeshConstraintsBlock
resource: sources/model_part_io.cpp#L5858-L5903
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/WriteInAllFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::DivideMeshConstraintsBlock( OutputFilesContainerType& OutputFiles, const PartitionIndicesContainerType& rConstraintsAllPartitions )`

# Calls

- [WriteInAllFiles](../../../../functions/sources/model_part_io/ModelPartIO/WriteInAllFiles.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReorderedConstraintId](../../../../functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId.md)

# Called by

- [DivideMeshBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshBlock.md)