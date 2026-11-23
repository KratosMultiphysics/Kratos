---
type: C++ Function
title: GetSolutionStepIndex
resource: includes/process_info.h#L277-L280
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/process_info/GetPreviousTimeStepIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/process_info/GetPreviousSolutionStepIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/process_info/ProcessInfo/CloneSolutionStepInfo
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/process_info/ProcessInfo/RemoveSolutionStepInfo
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`IndexType GetSolutionStepIndex() const`

# Called by

- [GetPreviousTimeStepIndex](../../../functions/includes/process_info/GetPreviousTimeStepIndex.md)
- [GetPreviousSolutionStepIndex](../../../functions/includes/process_info/GetPreviousSolutionStepIndex.md)
- [CloneSolutionStepInfo](../../../functions/sources/process_info/ProcessInfo/CloneSolutionStepInfo.md)
- [RemoveSolutionStepInfo](../../../functions/sources/process_info/ProcessInfo/RemoveSolutionStepInfo.md)