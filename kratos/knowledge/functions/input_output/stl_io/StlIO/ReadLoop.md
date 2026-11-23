---
type: C++ Method
title: ReadLoop
resource: input_output/stl_io.cpp#L436-L456
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/stl_io/StlIO/ReadKeyword
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/stl_io/StlIO/ReadFacet
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StlIO::ReadLoop( ModelPart & rThisModelPart, const std::function<void(ModelPart&, NodesArrayType&)>& rCreateEntityFunctor )`

# Calls

- [ReadKeyword](../../../../functions/input_output/stl_io/StlIO/ReadKeyword.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [ReadFacet](../../../../functions/input_output/stl_io/StlIO/ReadFacet.md)