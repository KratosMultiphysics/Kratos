---
type: Python Method
title: Save
resource: python_scripts/project.py#L69-L106
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def Save(self, save_folder_path: pathlib.Path, checkpoint_file_path: pathlib.Path, output_settings_file_path: Optional[pathlib.Path] = None) -> None:`

# Calls

- [MPISafeCreateDirectories](../../../../functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories.md)
- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)