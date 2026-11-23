---
type: C++ Method
title: ReadCommunicatorDataBlock
resource: sources/model_part_io.cpp#L3943-L3987
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/SetNumberOfColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorLocalNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorGhostNodesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadCommunicatorDataBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [SetNumberOfColors](../../../../functions/sources/communicator/Communicator/SetNumberOfColors.md)
- [ReadCommunicatorLocalNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorLocalNodesBlock.md)
- [ReadCommunicatorGhostNodesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorGhostNodesBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)

# Called by

- [ReadModelPart](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)