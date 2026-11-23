---
type: C++ Method
title: ReadNodesBlock
resource: sources/model_part_io.cpp#L1864-L1904
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/X0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Y0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Z0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadNodesBlock(NodesContainerType& rThisNodes)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [X0](../../../../functions/includes/node/X0.md)
- [Y0](../../../../functions/includes/node/Y0.md)
- [Z0](../../../../functions/includes/node/Z0.md)
- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)