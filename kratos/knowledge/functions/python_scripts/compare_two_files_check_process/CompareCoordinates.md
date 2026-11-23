---
type: Python Function
title: CompareCoordinates
resource: python_scripts/compare_two_files_check_process.py#L730-L767
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareGeoFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CompareCoordinates(lines_ref, lines_out, line_counter, error_messages):`

# Calls

- [__CheckCloseValues](../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CheckCloseValues.md)

# Called by

- [__CompareGeoFile](../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareGeoFile.md)