---
type: Python Function
title: GetNodalNonHistoricalVariableList
resource: mpi/python_scripts/debug_utilities.py#L64-L77
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/python_scripts/debug_utilities/GetVariableAndType
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/python_scripts/debug_utilities/GetNonHistoricalVariableList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetNodalNonHistoricalVariableList(modelPart, container, full_search=False, split_by_commponents=False):`

# Calls

- [GetVariableAndType](../../../../functions/mpi/python_scripts/debug_utilities/GetVariableAndType.md)

# Called by

- [GetNonHistoricalVariableList](../../../../functions/mpi/python_scripts/debug_utilities/GetNonHistoricalVariableList.md)