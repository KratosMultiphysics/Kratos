---
type: C++ Method
title: ReadNodalDofVariableData
resource: sources/model_part_io.cpp#L2833-L2866
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
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalDataBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadNodalDofVariableData(NodesContainerType& rThisNodes, const TVariableType& rVariable)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [FindKey](../../../../functions/sources/model_part_io/ModelPartIO/FindKey.md)
- [Fix](../../../../functions/includes/node/Fix.md)

# Called by

- [ReadNodalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalDataBlock.md)