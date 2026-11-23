---
type: C++ Method
title: SetCurrentTime
resource: sources/process_info.cpp#L63-L70
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/process_info/CreateTimeStepInfo
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/process_info/CloneTimeStepInfo
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/ReduceTimeStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/process_info/ProcessInfo/SetAsTimeStepInfo
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ProcessInfo::SetCurrentTime(double NewTime)`

# Called by

- [CreateTimeStepInfo](../../../../functions/includes/process_info/CreateTimeStepInfo.md)
- [CloneTimeStepInfo](../../../../functions/includes/process_info/CloneTimeStepInfo.md)
- [ReduceTimeStep](../../../../functions/sources/model_part/ModelPart/ReduceTimeStep.md)
- [SetAsTimeStepInfo](../../../../functions/sources/process_info/ProcessInfo/SetAsTimeStepInfo.md)