---
type: C++ Method
title: CloneTimeStep
resource: sources/model_part.cpp#L225-L235
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CloneSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPart::IndexType ModelPart::CloneTimeStep()`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [CloneSolutionStep](../../../../functions/sources/model_part/ModelPart/CloneSolutionStep.md)