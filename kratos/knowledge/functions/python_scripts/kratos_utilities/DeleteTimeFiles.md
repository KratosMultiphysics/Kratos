---
type: Python Function
title: DeleteTimeFiles
resource: python_scripts/kratos_utilities.py#L43-L46
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/DeleteFilesEndingWith
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/tearDown
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def DeleteTimeFiles(directory_name):`

# Calls

- [DeleteFilesEndingWith](../../../functions/python_scripts/kratos_utilities/DeleteFilesEndingWith.md)

# Called by

- [tearDown](../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/tearDown.md)