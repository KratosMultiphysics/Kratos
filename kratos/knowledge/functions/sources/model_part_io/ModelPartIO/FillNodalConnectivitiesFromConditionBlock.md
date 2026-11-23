---
type: C++ Method
title: FillNodalConnectivitiesFromConditionBlock
resource: sources/model_part_io.cpp#L3772-L3837
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalGraph
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::FillNodalConnectivitiesFromConditionBlock(ConnectivitiesContainerType& rNodalConnectivities)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadNodalGraph](../../../../functions/sources/model_part_io/ModelPartIO/ReadNodalGraph.md)