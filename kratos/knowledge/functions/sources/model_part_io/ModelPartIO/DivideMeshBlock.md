---
type: C++ Method
title: DivideMeshBlock
resource: sources/model_part_io.cpp#L5600-L5644
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
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshConstraintsBlock
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

`void ModelPartIO::DivideMeshBlock( OutputFilesContainerType& rOutputFiles, const PartitionIndicesContainerType& rNodesAllPartitions, const PartitionIndicesContainerType& rElementsAllPartitions, const PartitionIndicesContainerType& rConditionsAllPartitions, const PartitionIndicesContainerType& rConstraintsAllPartitions )`

# Calls

- [WriteInAllFiles](../../../../functions/sources/model_part_io/ModelPartIO/WriteInAllFiles.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [DivideMeshDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshDataBlock.md)
- [DivideMeshNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshNodesBlock.md)
- [DivideMeshElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshElementsBlock.md)
- [DivideMeshConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshConditionsBlock.md)
- [DivideMeshConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshConstraintsBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)

# Called by

- [DivideInputToPartitionsImpl](../../../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl.md)