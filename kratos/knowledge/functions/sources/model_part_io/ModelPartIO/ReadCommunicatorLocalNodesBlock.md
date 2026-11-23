---
type: C++ Method
title: ReadCommunicatorLocalNodesBlock
resource: sources/model_part_io.cpp#L3989-L4051
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/GetNumberOfColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FindKey
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadCommunicatorDataBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadCommunicatorLocalNodesBlock(Communicator& rThisCommunicator, NodesContainerType& rThisNodes)`

# Calls

- [GetNumberOfColors](../../../../functions/sources/communicator/Communicator/GetNumberOfColors.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [FindKey](../../../../functions/sources/model_part_io/ModelPartIO/FindKey.md)

# Called by

- [ReadCommunicatorDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadCommunicatorDataBlock.md)