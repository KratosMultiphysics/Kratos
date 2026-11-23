---
type: C++ Method
title: ReadVerticesAndFaces
resource: input_output/obj_io.cpp#L223-L269
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/obj_io/ObjIO/ParseVertexLine
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/obj_io/ObjIO/ParseNormalLine
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/obj_io/ObjIO/ParseFaceLine
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/obj_io/ObjIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ObjIO::ReadVerticesAndFaces( ModelPart& rThisModelPart, const std::string& rEntityType, const bool NormalAsHistoricalVariable, const bool DecomposeQuadrilateral )`

# Calls

- [ParseVertexLine](../../../../functions/input_output/obj_io/ObjIO/ParseVertexLine.md)
- [ParseNormalLine](../../../../functions/input_output/obj_io/ObjIO/ParseNormalLine.md)
- [ParseFaceLine](../../../../functions/input_output/obj_io/ObjIO/ParseFaceLine.md)

# Called by

- [ReadModelPart](../../../../functions/input_output/obj_io/ObjIO/ReadModelPart.md)