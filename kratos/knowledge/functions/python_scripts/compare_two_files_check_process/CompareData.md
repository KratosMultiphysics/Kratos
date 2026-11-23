---
type: Python Function
title: CompareData
resource: python_scripts/compare_two_files_check_process.py#L472-L526
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/compare_two_files_check_process/CompareFieldData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareVtkFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CompareData(lines_ref, lines_out, line_counter):`

# Calls

- [CompareFieldData](../../../functions/python_scripts/compare_two_files_check_process/CompareFieldData.md)

# Called by

- [__CompareVtkFile](../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareVtkFile.md)