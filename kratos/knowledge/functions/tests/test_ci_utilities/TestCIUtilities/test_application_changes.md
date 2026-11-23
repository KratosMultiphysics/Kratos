---
type: Python Method
title: test_application_changes
resource: tests/test_ci_utilities.py#L39-L47
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/testing/ci_utilities/get_changed_applications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/ci_utilities/is_core_changed
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/ci_utilities/is_mpi_core_changed
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/ci_utilities/get_changed_files_extensions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/ci_utilities/are_only_python_files_changed
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_application_changes(self, _):`

# Calls

- [get_changed_applications](../../../../functions/python_scripts/testing/ci_utilities/get_changed_applications.md)
- [is_core_changed](../../../../functions/python_scripts/testing/ci_utilities/is_core_changed.md)
- [is_mpi_core_changed](../../../../functions/python_scripts/testing/ci_utilities/is_mpi_core_changed.md)
- [get_changed_files_extensions](../../../../functions/python_scripts/testing/ci_utilities/get_changed_files_extensions.md)
- [are_only_python_files_changed](../../../../functions/python_scripts/testing/ci_utilities/are_only_python_files_changed.md)