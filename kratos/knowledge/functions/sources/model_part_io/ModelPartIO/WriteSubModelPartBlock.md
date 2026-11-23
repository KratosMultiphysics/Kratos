---
type: C++ Method
title: WriteSubModelPartBlock
resource: sources/model_part_io.cpp#L4485-L4558
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/GetSubModelPartNames
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/WriteModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::WriteSubModelPartBlock( ModelPart& rMainModelPart, const std::string& InitialTabulation)`

# Calls

- [GetSubModelPartNames](../../../../functions/sources/model_part/ModelPart/GetSubModelPartNames.md)

# Called by

- [WriteModelPart](../../../../functions/sources/model_part_io/ModelPartIO/WriteModelPart.md)