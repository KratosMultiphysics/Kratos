---
type: C++ Method
title: CreateTimeStep
resource: sources/model_part.cpp#L238-L248
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPart::IndexType ModelPart::CreateTimeStep(double NewTime)`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [CreateSolutionStep](../../../../functions/sources/model_part/ModelPart/CreateSolutionStep.md)