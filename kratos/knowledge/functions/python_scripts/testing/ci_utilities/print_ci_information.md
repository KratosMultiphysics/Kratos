---
type: Python Function
title: print_ci_information
resource: python_scripts/testing/ci_utilities.py#L63-L71
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/testing/ci_utilities/changed_files
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/ci_utilities/ci_applications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/ci_utilities/get_changed_files_extensions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/ci_utilities/are_only_python_files_changed
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/ci_utilities/get_changed_applications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/ci_utilities/is_core_changed
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/ci_utilities/is_mpi_core_changed
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def print_ci_information() -> None:`

# Calls

- [changed_files](../../../../functions/python_scripts/testing/ci_utilities/changed_files.md)
- [ci_applications](../../../../functions/python_scripts/testing/ci_utilities/ci_applications.md)
- [get_changed_files_extensions](../../../../functions/python_scripts/testing/ci_utilities/get_changed_files_extensions.md)
- [are_only_python_files_changed](../../../../functions/python_scripts/testing/ci_utilities/are_only_python_files_changed.md)
- [get_changed_applications](../../../../functions/python_scripts/testing/ci_utilities/get_changed_applications.md)
- [is_core_changed](../../../../functions/python_scripts/testing/ci_utilities/is_core_changed.md)
- [is_mpi_core_changed](../../../../functions/python_scripts/testing/ci_utilities/is_mpi_core_changed.md)