---
type: C++ Method
title: ReadNodalGraph
resource: sources/model_part_io.cpp#L945-L1012
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ResetInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ScanNodeBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromGeometryBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t ModelPartIO::ReadNodalGraph(ConnectivitiesContainerType& rAuxConnectivities)`

# Calls

- [ResetInput](../../../../functions/sources/model_part_io/ModelPartIO/ResetInput.md)
- [ScanNodeBlock](../../../../functions/sources/model_part_io/ModelPartIO/ScanNodeBlock.md)
- [FillNodalConnectivitiesFromGeometryBlock](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromGeometryBlock.md)
- [FillNodalConnectivitiesFromElementBlock](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlock.md)
- [FillNodalConnectivitiesFromConditionBlock](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlock.md)
- [FillNodalConnectivitiesFromConstraintBlock](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)