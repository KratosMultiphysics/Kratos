---
type: C++ Function
title: CloneTimeStepInfo
resource: includes/process_info.h#L131-L135
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/process_info/ProcessInfo/CreateSolutionStepInfo
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/process_info/ProcessInfo/SetCurrentTime
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CloneTimeStepInfo(IndexType SolutionStepIndex = 0)`

# Calls

- [CreateSolutionStepInfo](../../../functions/sources/process_info/ProcessInfo/CreateSolutionStepInfo.md)
- [SetCurrentTime](../../../functions/sources/process_info/ProcessInfo/SetCurrentTime.md)