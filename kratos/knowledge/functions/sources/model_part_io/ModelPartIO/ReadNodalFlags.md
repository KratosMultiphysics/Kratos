---
type: C++ Method
title: ReadNodalFlags
resource: sources/model_part_io.cpp#L2868-L2889
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FindKey
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalDataBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadNodalFlags(NodesContainerType& rThisNodes, Flags const& rFlags)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [FindKey](../../../../functions/sources/model_part_io/ModelPartIO/FindKey.md)

# Called by

- [ReadNodalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalDataBlock.md)