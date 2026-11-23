---
type: C++ Method
title: CreateSolutionStepInfo
resource: sources/process_info.cpp#L80-L88
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
  - target: functions/sources/process_info/ProcessInfo/CloneSolutionStepInfo
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ProcessInfo::CreateSolutionStepInfo(IndexType SolutionStepIndex)`

# Called by

- [CreateTimeStepInfo](../../../../functions/includes/process_info/CreateTimeStepInfo.md)
- [CloneTimeStepInfo](../../../../functions/includes/process_info/CloneTimeStepInfo.md)
- [CloneSolutionStepInfo](../../../../functions/sources/process_info/ProcessInfo/CloneSolutionStepInfo.md)