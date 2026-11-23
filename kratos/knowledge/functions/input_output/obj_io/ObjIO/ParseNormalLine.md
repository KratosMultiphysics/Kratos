---
type: C++ Method
title: ParseNormalLine
resource: input_output/obj_io.cpp#L295-L321
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

`void ObjIO::ParseNormalLine( ModelPart& rThisModelPart, const std::string& rLine, const bool NormalAsHistoricalVariable )`

# Calls

- [Tokenize](../../../../functions/input_output/obj_io/ObjIO/Tokenize.md)

# Called by

- [ReadVerticesAndFaces](../../../../functions/input_output/obj_io/ObjIO/ReadVerticesAndFaces.md)