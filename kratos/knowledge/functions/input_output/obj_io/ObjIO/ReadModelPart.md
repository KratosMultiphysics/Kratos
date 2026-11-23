---
type: C++ Method
title: ReadModelPart
resource: input_output/obj_io.cpp#L97-L125
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
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/obj_io/ObjIO/ReadVerticesAndFaces
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ObjIO::ReadModelPart(ModelPart& rThisModelPart)`

# Calls

- [RecursivelyHasProperties](../../../../functions/sources/model_part/ModelPart/RecursivelyHasProperties.md)
- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ReadVerticesAndFaces](../../../../functions/input_output/obj_io/ObjIO/ReadVerticesAndFaces.md)