---
type: Python Function
title: are_only_python_files_changed
resource: python_scripts/testing/ci_utilities.py#L59-L60
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/testing/ci_utilities/get_changed_files_extensions
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/testing/ci_utilities/print_ci_information
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ci_utilities/TestCIUtilities/test_core_changes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ci_utilities/TestCIUtilities/test_mpi_core_changes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ci_utilities/TestCIUtilities/test_application_changes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def are_only_python_files_changed() -> bool:`

# Calls

- [get_changed_files_extensions](../../../../functions/python_scripts/testing/ci_utilities/get_changed_files_extensions.md)

# Called by

- [print_ci_information](../../../../functions/python_scripts/testing/ci_utilities/print_ci_information.md)
- [test_core_changes](../../../../functions/tests/test_ci_utilities/TestCIUtilities/test_core_changes.md)
- [test_mpi_core_changes](../../../../functions/tests/test_ci_utilities/TestCIUtilities/test_mpi_core_changes.md)
- [test_application_changes](../../../../functions/tests/test_ci_utilities/TestCIUtilities/test_application_changes.md)