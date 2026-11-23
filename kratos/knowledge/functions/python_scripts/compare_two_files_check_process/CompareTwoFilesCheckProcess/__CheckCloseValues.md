---
type: Python Method
title: __CheckCloseValues
resource: python_scripts/compare_two_files_check_process.py#L968-L975
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__ComparePostResFileResultsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDelimittedFileResults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDatFileResultsWithLocation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/ComparePoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareFieldData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareCaseFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareCoordinates
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CheckCloseValues(self, val_a, val_b, additional_info=""):`

# Called by

- [__ComparePostResFileResultsBlock](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__ComparePostResFileResultsBlock.md)
- [__CompareDelimittedFileResults](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDelimittedFileResults.md)
- [__CompareDatFileResultsWithLocation](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareDatFileResultsWithLocation.md)
- [ComparePoints](../../../../functions/python_scripts/compare_two_files_check_process/ComparePoints.md)
- [CompareFieldData](../../../../functions/python_scripts/compare_two_files_check_process/CompareFieldData.md)
- [__CompareCaseFile](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__CompareCaseFile.md)
- [CompareCoordinates](../../../../functions/python_scripts/compare_two_files_check_process/CompareCoordinates.md)