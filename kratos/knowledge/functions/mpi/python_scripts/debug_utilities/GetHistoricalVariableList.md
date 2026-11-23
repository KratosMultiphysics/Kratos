---
type: Python Function
title: GetHistoricalVariableList
resource: mpi/python_scripts/debug_utilities.py#L60-L62
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/python_scripts/debug_utilities/GetVariableAndType
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/python_scripts/debug_utilities/CheckAllHistoricalVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/python_scripts/debug_utilities/CheckAllNonHistoricalVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/python_scripts/debug_utilities/CheckAllHistoricalVariablesStatistics
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/python_scripts/debug_utilities/CheckAllNonHistoricalVariablesStatistics
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetHistoricalVariableList(model_part, split_by_commponents=False):`

# Calls

- [GetVariableAndType](../../../../functions/mpi/python_scripts/debug_utilities/GetVariableAndType.md)

# Called by

- [CheckAllHistoricalVariables](../../../../functions/mpi/python_scripts/debug_utilities/CheckAllHistoricalVariables.md)
- [CheckAllNonHistoricalVariables](../../../../functions/mpi/python_scripts/debug_utilities/CheckAllNonHistoricalVariables.md)
- [CheckAllHistoricalVariablesStatistics](../../../../functions/mpi/python_scripts/debug_utilities/CheckAllHistoricalVariablesStatistics.md)
- [CheckAllNonHistoricalVariablesStatistics](../../../../functions/mpi/python_scripts/debug_utilities/CheckAllNonHistoricalVariablesStatistics.md)