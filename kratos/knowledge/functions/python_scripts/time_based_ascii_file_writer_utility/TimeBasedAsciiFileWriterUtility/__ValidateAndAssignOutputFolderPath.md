---
type: Python Method
title: __ValidateAndAssignOutputFolderPath
resource: python_scripts/time_based_ascii_file_writer_utility.py#L131-L148
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __ValidateAndAssignOutputFolderPath(self): # a file name must be specified`

# Calls

- [MPISafeCreateDirectories](../../../../functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories.md)

# Called by

- [__init__](../../../../functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__init__.md)