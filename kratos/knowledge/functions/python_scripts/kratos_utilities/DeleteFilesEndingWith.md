---
type: Python Function
title: DeleteFilesEndingWith
resource: python_scripts/kratos_utilities.py#L34-L41
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/DeleteFileIfExisting
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/kratos_utilities/DeleteTimeFiles
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def DeleteFilesEndingWith(directory_name, *ends_with):`

# Calls

- [DeleteFileIfExisting](../../../functions/python_scripts/kratos_utilities/DeleteFileIfExisting.md)

# Called by

- [DeleteTimeFiles](../../../functions/python_scripts/kratos_utilities/DeleteTimeFiles.md)