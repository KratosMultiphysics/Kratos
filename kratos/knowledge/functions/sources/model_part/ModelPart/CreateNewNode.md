---
type: C++ Method
title: CreateNewNode
resource: sources/model_part.cpp#L324-L364
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/SetSolutionStepVariablesList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPart::NodeType::Pointer ModelPart::CreateNewNode(int Id, double x, double y, double z, VariablesList::Pointer pNewVariablesList, ModelPart::IndexType ThisIndex)`

# Calls

- [SetSolutionStepVariablesList](../../../../functions/includes/node/SetSolutionStepVariablesList.md)