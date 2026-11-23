---
type: C++ Method
title: Tokenize
resource: input_output/obj_io.cpp#L394-L419
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/pointer_vector/emplace_back
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/obj_io/ObjIO/ParseVertexLine
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/obj_io/ObjIO/ParseNormalLine
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/obj_io/ObjIO/ParseFaceLine
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::string> ObjIO::Tokenize(const std::string& rLine)`

# Calls

- [emplace_back](../../../../functions/containers/pointer_vector/emplace_back.md)

# Called by

- [ParseVertexLine](../../../../functions/input_output/obj_io/ObjIO/ParseVertexLine.md)
- [ParseNormalLine](../../../../functions/input_output/obj_io/ObjIO/ParseNormalLine.md)
- [ParseFaceLine](../../../../functions/input_output/obj_io/ObjIO/ParseFaceLine.md)