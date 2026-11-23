---
type: C++ Method
title: CloneSolutionStep
resource: sources/model_part.cpp#L206-L223
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/CloneSolutionStepData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/process_info/ProcessInfo/ClearHistory
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part/ModelPart/CloneTimeStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPart::IndexType ModelPart::CloneSolutionStep()`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [CloneSolutionStepData](../../../../functions/includes/node/CloneSolutionStepData.md)
- [ClearHistory](../../../../functions/sources/process_info/ProcessInfo/ClearHistory.md)

# Called by

- [CloneTimeStep](../../../../functions/sources/model_part/ModelPart/CloneTimeStep.md)