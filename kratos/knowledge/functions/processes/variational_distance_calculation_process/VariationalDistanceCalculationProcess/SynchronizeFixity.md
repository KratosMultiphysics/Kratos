---
type: C++ Method
title: SynchronizeFixity
resource: processes/variational_distance_calculation_process.h#L675-L696
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/TotalProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SynchronizeFixity()`

# Calls

- [TotalProcesses](../../../../functions/sources/communicator/Communicator/TotalProcesses.md)
- [Fix](../../../../functions/includes/node/Fix.md)

# Called by

- [Execute](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute.md)