---
type: Python Method
title: __AddToExistingOutputFile
resource: python_scripts/time_based_ascii_file_writer_utility.py#L79-L129
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__OpenOutputFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __AddToExistingOutputFile(self, file_header, restart_time):`

# Calls

- [__OpenOutputFile](../../../../functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__OpenOutputFile.md)

# Called by

- [__init__](../../../../functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__init__.md)