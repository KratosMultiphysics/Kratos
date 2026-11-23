---
type: C++ Method
title: WriteEntityBlock
resource: input_output/stl_io.cpp#L178-L194
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/IsValidGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/WriteFacet
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/stl_io/StlIO/WriteModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StlIO::WriteEntityBlock(const TContainerType& rThisEntities)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [IsValidGeometry](../../../../functions/input_output/stl_io/StlIO/IsValidGeometry.md)
- [WriteFacet](../../../../functions/input_output/stl_io/StlIO/WriteFacet.md)

# Called by

- [WriteModelPart](../../../../functions/input_output/stl_io/StlIO/WriteModelPart.md)