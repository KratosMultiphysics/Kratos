---
type: C++ Method
title: ResetInput
resource: sources/model_part_io.cpp#L6557-L6562
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodesNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionsConnectivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadInitialValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalGraph
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ResetInput()`

# Called by

- [ReadNodesNumber](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodesNumber.md)
- [ReadGeometriesConnectivities](../../../../functions/sources/model_part_io/ModelPartIO/ReadGeometriesConnectivities.md)
- [ReadElementsConnectivities](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivities.md)
- [ReadConditionsConnectivities](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionsConnectivities.md)
- [ReadConstraintsConnectivities](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsConnectivities.md)
- [ReadInitialValues](../../../../functions/sources/model_part_io/ModelPartIO/ReadInitialValues.md)
- [ReadModelPart](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)
- [ReadNodalGraph](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraph.md)