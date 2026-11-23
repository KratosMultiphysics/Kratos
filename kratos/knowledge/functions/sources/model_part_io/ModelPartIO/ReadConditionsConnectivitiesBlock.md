---
type: C++ Method
title: ReadConditionsConnectivitiesBlock
resource: sources/model_part_io.cpp#L3470-L3530
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
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionsConnectivities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPartIO::SizeType ModelPartIO::ReadConditionsConnectivitiesBlock(ConnectivitiesContainerType& rThisConnectivities)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadConditionsConnectivities](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionsConnectivities.md)