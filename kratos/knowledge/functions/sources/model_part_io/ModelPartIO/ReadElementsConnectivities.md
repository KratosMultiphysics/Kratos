---
type: C++ Method
title: ReadElementsConnectivities
resource: sources/model_part_io.cpp#L354-L374
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ResetInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivitiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t ModelPartIO::ReadElementsConnectivities(ConnectivitiesContainerType& rElementsConnectivities)`

# Calls

- [ResetInput](../../../../functions/sources/model_part_io/ModelPartIO/ResetInput.md)
- [ReadElementsConnectivitiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivitiesBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)