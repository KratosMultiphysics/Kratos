---
type: C++ Method
title: DivideSubModelPartBlock
resource: sources/model_part_io.cpp#L5646-L5695
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
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartTableBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartGeometriesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::DivideSubModelPartBlock( OutputFilesContainerType& rOutputFiles, const PartitionIndicesContainerType& rNodesAllPartitions, const PartitionIndicesContainerType& rElementsAllPartitions, const PartitionIndicesContainerType& rConditionsAllPartitions, const PartitionIndicesContainerType& rConstraintsAllPartitions, const PartitionIndicesContainerType& rGeometriesAllPartitions )`

# Calls

- [WriteInAllFiles](../../../../functions/sources/model_part_io/ModelPartIO/WriteInAllFiles.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [DivideSubModelPartDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartDataBlock.md)
- [DivideSubModelPartTableBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartTableBlock.md)
- [DivideSubModelPartNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartNodesBlock.md)
- [DivideSubModelPartElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartElementsBlock.md)
- [DivideSubModelPartConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartConditionsBlock.md)
- [DivideSubModelPartConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartConstraintsBlock.md)
- [DivideSubModelPartGeometriesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartGeometriesBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)

# Called by

- [DivideInputToPartitionsImpl](../../../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl.md)