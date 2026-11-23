---
type: C++ Method
title: CalculateTimeIntegratedNonHistoricalQuantity
resource: processes/time_averaging_process.cpp#L385-L402
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/time_averaging_process/if
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TimeAveragingProcess::CalculateTimeIntegratedNonHistoricalQuantity( TContainerType& rContainer, const Variable<TDataType>& rVariable, const Variable<TDataType>& rAveragedVariable, const double DeltaTime) const`

# Called by

- [if](../../../../functions/processes/time_averaging_process/if.md)