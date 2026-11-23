---
type: C++ Method
title: ReadConstraintsConnectivities
resource: sources/model_part_io.cpp#L529-L547
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ResetInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivitiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t ModelPartIO::ReadConstraintsConnectivities(ConnectivitiesContainerType& rConstraintsConnectivities)`

# Calls

- [ResetInput](../../../../functions/sources/model_part_io/ModelPartIO/ResetInput.md)
- [ReadConstraintsConnectivitiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivitiesBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)