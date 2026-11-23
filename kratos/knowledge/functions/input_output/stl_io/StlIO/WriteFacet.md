---
type: C++ Method
title: WriteFacet
resource: input_output/stl_io.cpp#L320-L336
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/stl_io/StlIO/WriteEntityBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/WriteGeometryBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/WriteEntityBlockMPI
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/WriteGeometryBlockMPI
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StlIO::WriteFacet( const GeometryType& rGeom, TStreamType& rStream )`

# Called by

- [WriteEntityBlock](../../../../functions/input_output/stl_io/StlIO/WriteEntityBlock.md)
- [WriteGeometryBlock](../../../../functions/input_output/stl_io/StlIO/WriteGeometryBlock.md)
- [WriteEntityBlockMPI](../../../../functions/input_output/stl_io/StlIO/WriteEntityBlockMPI.md)
- [WriteGeometryBlockMPI](../../../../functions/input_output/stl_io/StlIO/WriteGeometryBlockMPI.md)