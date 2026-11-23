---
type: C++ Method
title: ReadElementsConnectivitiesBlock
resource: sources/model_part_io.cpp#L3410-L3468
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
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPartIO::SizeType ModelPartIO::ReadElementsConnectivitiesBlock(ConnectivitiesContainerType& rThisConnectivities)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadElementsConnectivities](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivities.md)