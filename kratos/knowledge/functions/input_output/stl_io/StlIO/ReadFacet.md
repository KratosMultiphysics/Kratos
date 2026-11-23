---
type: C++ Method
title: ReadFacet
resource: input_output/stl_io.cpp#L415-L434
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/stl_io/StlIO/ReadKeyword
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/ReadLoop
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/stl_io/StlIO/ReadSolid
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StlIO::ReadFacet( ModelPart& rThisModelPart, const std::function<void(ModelPart&, NodesArrayType&)>& rCreateEntityFunctor )`

# Calls

- [ReadKeyword](../../../../functions/input_output/stl_io/StlIO/ReadKeyword.md)
- [ReadLoop](../../../../functions/input_output/stl_io/StlIO/ReadLoop.md)

# Called by

- [ReadSolid](../../../../functions/input_output/stl_io/StlIO/ReadSolid.md)