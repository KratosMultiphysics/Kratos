---
type: C++ Function
title: CreateTimeStepInfo
resource: includes/process_info.h#L125-L129
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

`void CreateTimeStepInfo(IndexType SolutionStepIndex = 0)`

# Calls

- [CreateSolutionStepInfo](../../../functions/sources/process_info/ProcessInfo/CreateSolutionStepInfo.md)
- [SetCurrentTime](../../../functions/sources/process_info/ProcessInfo/SetCurrentTime.md)