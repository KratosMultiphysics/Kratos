---
type: C++ Method
title: ReadNodalScalarVariableData
resource: sources/model_part_io.cpp#L2892-L2929
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

`void ModelPartIO::ReadNodalScalarVariableData(NodesContainerType& rThisNodes, const TVariableType& rVariable)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [FindKey](../../../../functions/sources/model_part_io/ModelPartIO/FindKey.md)

# Called by

- [ReadNodalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalDataBlock.md)