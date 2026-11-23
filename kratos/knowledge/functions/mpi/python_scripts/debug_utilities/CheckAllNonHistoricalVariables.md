---
type: Python Function
title: CheckAllNonHistoricalVariables
resource: mpi/python_scripts/debug_utilities.py#L107-L116
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/python_scripts/debug_utilities/GetHistoricalVariableList
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/python_scripts/debug_utilities/CheckAll
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CheckAllNonHistoricalVariables(model_part):`

# Calls

- [GetHistoricalVariableList](../../../../functions/mpi/python_scripts/debug_utilities/GetHistoricalVariableList.md)

# Called by

- [CheckAll](../../../../functions/mpi/python_scripts/debug_utilities/CheckAll.md)