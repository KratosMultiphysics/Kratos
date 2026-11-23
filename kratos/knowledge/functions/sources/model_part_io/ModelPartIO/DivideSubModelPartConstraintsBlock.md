---
type: C++ Method
title: DivideSubModelPartConstraintsBlock
resource: sources/model_part_io.cpp#L6080-L6126
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
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::DivideSubModelPartConstraintsBlock( OutputFilesContainerType& rOutputFiles, const PartitionIndicesContainerType& rConstraintsAllPartitions )`

# Calls

- [WriteInAllFiles](../../../../functions/sources/model_part_io/ModelPartIO/WriteInAllFiles.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReorderedConstraintId](../../../../functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId.md)

# Called by

- [DivideSubModelPartBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartBlock.md)