---
type: C++ Function
title: if
resource: processes/time_averaging_process.cpp#L196-L269
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/time_averaging_process/TimeAveragingProcess/CalculateTimeIntegratedNodalHistoricalQuantity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/time_averaging_process/TimeAveragingProcess/CalculateTimeIntegratedNonHistoricalQuantity
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TRY if (IsIntegrationStep())`

# Calls

- [CalculateTimeIntegratedNodalHistoricalQuantity](../../../functions/processes/time_averaging_process/TimeAveragingProcess/CalculateTimeIntegratedNodalHistoricalQuantity.md)
- [CalculateTimeIntegratedNonHistoricalQuantity](../../../functions/processes/time_averaging_process/TimeAveragingProcess/CalculateTimeIntegratedNonHistoricalQuantity.md)