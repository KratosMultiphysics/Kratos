---
type: C++ Method
title: WriteModelPart
resource: input_output/stl_io.cpp#L133-L175
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/stl_io/StlIO/WriteEntityBlockMPI
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/WriteGeometryBlockMPI
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/WriteEntityBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/WriteGeometryBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StlIO::WriteModelPart(const ModelPart& rThisModelPart)`

# Calls

- [WriteEntityBlockMPI](../../../../functions/input_output/stl_io/StlIO/WriteEntityBlockMPI.md)
- [WriteGeometryBlockMPI](../../../../functions/input_output/stl_io/StlIO/WriteGeometryBlockMPI.md)
- [WriteEntityBlock](../../../../functions/input_output/stl_io/StlIO/WriteEntityBlock.md)
- [WriteGeometryBlock](../../../../functions/input_output/stl_io/StlIO/WriteGeometryBlock.md)