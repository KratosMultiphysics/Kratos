---
type: C++ Method
title: ParseFaceLine
resource: input_output/obj_io.cpp#L326-L389
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/obj_io/ObjIO/Tokenize
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/obj_io/ObjIO/ReadVerticesAndFaces
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ObjIO::ParseFaceLine( ModelPart& rThisModelPart, const std::string& rLine, const std::string& rEntityType, const bool DecomposeQuadrilateral )`

# Calls

- [Tokenize](../../../../functions/input_output/obj_io/ObjIO/Tokenize.md)

# Called by

- [ReadVerticesAndFaces](../../../../functions/input_output/obj_io/ObjIO/ReadVerticesAndFaces.md)