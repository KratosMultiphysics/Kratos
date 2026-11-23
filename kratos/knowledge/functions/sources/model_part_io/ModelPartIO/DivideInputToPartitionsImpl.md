---
type: C++ Method
title: DivideInputToPartitionsImpl
resource: sources/model_part_io.cpp#L1479-L1548
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/DivideModelPartDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideTableBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DividePropertiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideGeometriesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideNodalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideElementalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideConditionalDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideConstraintDataBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideMeshBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideSubModelPartBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WritePartitionIndices
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteCommunicatorData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/DivideInputToPartitions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::DivideInputToPartitionsImpl( OutputFilesContainerType& rOutputFiles, SizeType NumberOfPartitions, const PartitioningInfo& rPartitioningInfo )`

# Calls

- [DivideModelPartDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideModelPartDataBlock.md)
- [DivideTableBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideTableBlock.md)
- [DividePropertiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DividePropertiesBlock.md)
- [DivideNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideNodesBlock.md)
- [DivideGeometriesBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideGeometriesBlock.md)
- [DivideElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideElementsBlock.md)
- [DivideConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideConditionsBlock.md)
- [DivideConstraintsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideConstraintsBlock.md)
- [DivideNodalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideNodalDataBlock.md)
- [DivideElementalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideElementalDataBlock.md)
- [DivideConditionalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideConditionalDataBlock.md)
- [DivideConstraintDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideConstraintDataBlock.md)
- [DivideMeshBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideMeshBlock.md)
- [DivideSubModelPartBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideSubModelPartBlock.md)
- [WritePartitionIndices](../../../../functions/sources/model_part_io/ModelPartIO/WritePartitionIndices.md)
- [WriteCommunicatorData](../../../../functions/sources/model_part_io/ModelPartIO/WriteCommunicatorData.md)

# Called by

- [DivideInputToPartitions](../../../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitions.md)