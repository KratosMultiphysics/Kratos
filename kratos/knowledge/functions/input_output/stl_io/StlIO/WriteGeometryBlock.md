---
type: C++ Method
title: WriteGeometryBlock
resource: input_output/stl_io.cpp#L196-L211
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
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

`void StlIO::WriteGeometryBlock(const GeometryContainerType& rThisGeometries)`

# Calls

- [IsValidGeometry](../../../../functions/input_output/stl_io/StlIO/IsValidGeometry.md)
- [WriteFacet](../../../../functions/input_output/stl_io/StlIO/WriteFacet.md)

# Called by

- [WriteModelPart](../../../../functions/input_output/stl_io/StlIO/WriteModelPart.md)