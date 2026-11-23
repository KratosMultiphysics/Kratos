---
type: C++ Method
title: ReadGeometriesConnectivities
resource: sources/model_part_io.cpp#L256-L275
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ResetInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivitiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t ModelPartIO::ReadGeometriesConnectivities(ConnectivitiesContainerType& rGeometriesConnectivities)`

# Calls

- [ResetInput](../../../../functions/sources/model_part_io/ModelPartIO/ResetInput.md)
- [ReadGeometriesConnectivitiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivitiesBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)