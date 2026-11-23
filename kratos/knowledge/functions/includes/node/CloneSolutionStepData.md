---
type: C++ Function
title: CloneSolutionStepData
resource: includes/node.h#L494-L497
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/variables_list_data_value_container/CloneFront
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part/ModelPart/CloneSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CloneSolutionStepData()`

# Calls

- [CloneFront](../../../functions/containers/variables_list_data_value_container/CloneFront.md)

# Called by

- [CloneSolutionStep](../../../functions/sources/model_part/ModelPart/CloneSolutionStep.md)