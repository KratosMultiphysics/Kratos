---
type: C++ Method
title: ModelPart
resource: sources/model_part.cpp#L118
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/SetLocalMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPart::ModelPart(VariablesList::Pointer pVariablesList, Model& rOwnerModel) : ModelPart("Default", pVariablesList, rOwnerModel)`

# Calls

- [SetLocalMesh](../../../../functions/sources/communicator/Communicator/SetLocalMesh.md)