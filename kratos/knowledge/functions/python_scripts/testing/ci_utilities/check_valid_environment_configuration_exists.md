---
type: Python Function
title: check_valid_environment_configuration_exists
resource: python_scripts/testing/ci_utilities.py#L9-L22
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_filesystem/exists
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/testing/ci_utilities/changed_files
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/ci_utilities/ci_applications
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def check_valid_environment_configuration_exists() -> None:`

# Calls

- [exists](../../../../functions/sources/kratos_filesystem/exists.md)

# Called by

- [changed_files](../../../../functions/python_scripts/testing/ci_utilities/changed_files.md)
- [ci_applications](../../../../functions/python_scripts/testing/ci_utilities/ci_applications.md)