---
type: C++ Method
title: ReadModelPart
resource: input_output/stl_io.cpp#L86-L131
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/RecursivelyHasProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/ReadSolid
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StlIO::ReadModelPart(ModelPart& rThisModelPart)`

# Calls

- [RecursivelyHasProperties](../../../../functions/sources/model_part/ModelPart/RecursivelyHasProperties.md)
- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ReadSolid](../../../../functions/input_output/stl_io/StlIO/ReadSolid.md)