---
type: C++ Method
title: RootMeanSquareMethod
resource: processes/time_averaging_process.cpp#L429-L436
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

`void TimeAveragingProcess::RootMeanSquareMethod(const TDataType& rTemporalVariable, TDataType& rAveragedVariable, const double DeltaTime) const`

# Called by

- [GetTimeAveragingMethod](../../../../functions/processes/time_averaging_process/TimeAveragingProcess/GetTimeAveragingMethod.md)