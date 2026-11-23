---
type: C++ Method
title: GetTimeAveragingMethod
resource: processes/time_averaging_process.cpp#L341-L365
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/time_averaging_process/TimeAveragingProcess/AverageMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/time_averaging_process/TimeAveragingProcess/RootMeanSquareMethod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::function<void(const TDataType&, TDataType&, const double)> TimeAveragingProcess::GetTimeAveragingMethod() const`

# Calls

- [AverageMethod](../../../../functions/processes/time_averaging_process/TimeAveragingProcess/AverageMethod.md)
- [RootMeanSquareMethod](../../../../functions/processes/time_averaging_process/TimeAveragingProcess/RootMeanSquareMethod.md)