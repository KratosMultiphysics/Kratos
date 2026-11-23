---
type: C++ Method
title: AverageMethod
resource: processes/time_averaging_process.cpp#L405-L411
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/time_averaging_process/TimeAveragingProcess/GetTimeAveragingMethod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TimeAveragingProcess::AverageMethod(const TDataType& rTemporalVariable, TDataType& rAveragedVariable, const double DeltaTime) const`

# Called by

- [GetTimeAveragingMethod](../../../../functions/processes/time_averaging_process/TimeAveragingProcess/GetTimeAveragingMethod.md)