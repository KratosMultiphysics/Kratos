---
type: C++ Method
title: IsValidGeometry
resource: input_output/stl_io.cpp#L354-L368
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

`bool StlIO::IsValidGeometry( const Geometry<Node>& rGeometry, IndexType& rNumDegenerateGeos ) const`

# Called by

- [WriteEntityBlock](../../../../functions/input_output/stl_io/StlIO/WriteEntityBlock.md)
- [WriteGeometryBlock](../../../../functions/input_output/stl_io/StlIO/WriteGeometryBlock.md)
- [WriteEntityBlockMPI](../../../../functions/input_output/stl_io/StlIO/WriteEntityBlockMPI.md)
- [WriteGeometryBlockMPI](../../../../functions/input_output/stl_io/StlIO/WriteGeometryBlockMPI.md)