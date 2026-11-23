---
type: Python Function
title: CompareFieldData
resource: python_scripts/compare_two_files_check_process.py#L473-L501
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CompareFieldData(lines_ref, lines_out, line_counter, num_entities):`

# Calls

- [__CheckCloseValues](../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues.md)

# Called by

- [CompareData](../../../functions/python_scripts/compare_two_files_check_process/CompareData.md)