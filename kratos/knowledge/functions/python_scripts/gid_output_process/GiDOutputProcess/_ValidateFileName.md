---
type: Python Method
title: _ValidateFileName
resource: python_scripts/gid_output_process.py#L737-L746
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _ValidateFileName(file_name): # A file name must be specified`

# Calls

- [MPISafeCreateDirectories](../../../../functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories.md)

# Called by

- [__init__](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__init__.md)