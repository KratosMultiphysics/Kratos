---
type: C++ Method
title: CloneSolutionStepInfo
resource: sources/process_info.cpp#L90-L97
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/process_info/GetSolutionStepIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/process_info/ProcessInfo/CreateSolutionStepInfo
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ProcessInfo::CloneSolutionStepInfo()`

# Calls

- [GetSolutionStepIndex](../../../../functions/includes/process_info/GetSolutionStepIndex.md)
- [CreateSolutionStepInfo](../../../../functions/sources/process_info/ProcessInfo/CreateSolutionStepInfo.md)