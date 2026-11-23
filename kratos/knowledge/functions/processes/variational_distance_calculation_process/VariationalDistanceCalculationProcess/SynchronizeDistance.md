---
type: C++ Method
title: SynchronizeDistance
resource: processes/variational_distance_calculation_process.h#L646-L673
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/TotalProcesses
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SynchronizeDistance()`

# Calls

- [TotalProcesses](../../../../functions/sources/communicator/Communicator/TotalProcesses.md)

# Called by

- [Execute](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute.md)