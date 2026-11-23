---
type: Python Method
title: __init__
resource: python_scripts/time_based_ascii_file_writer_utility.py#L13-L69
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RemoveValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__ValidateAndAssignOutputFolderPath
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/TotalProcesses
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__InitializeOutputFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__AddToExistingOutputFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, model_part, params, file_header):`

# Calls

- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)
- [RemoveValue](../../../../functions/sources/kratos_parameters/Parameters/RemoveValue.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [__ValidateAndAssignOutputFolderPath](../../../../functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__ValidateAndAssignOutputFolderPath.md)
- [TotalProcesses](../../../../functions/sources/communicator/Communicator/TotalProcesses.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [__InitializeOutputFile](../../../../functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__InitializeOutputFile.md)
- [__AddToExistingOutputFile](../../../../functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__AddToExistingOutputFile.md)