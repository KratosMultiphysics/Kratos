---
type: Python Function
title: is_mpi_core_changed
resource: python_scripts/testing/ci_utilities.py#L51-L52
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/testing/ci_utilities/changed_files
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

`def is_mpi_core_changed() -> bool:`

# Calls

- [changed_files](../../../../functions/python_scripts/testing/ci_utilities/changed_files.md)

# Called by

- [print_ci_information](../../../../functions/python_scripts/testing/ci_utilities/print_ci_information.md)
- [test_core_changes](../../../../functions/tests/test_ci_utilities/TestCIUtilities/test_core_changes.md)
- [test_mpi_core_changes](../../../../functions/tests/test_ci_utilities/TestCIUtilities/test_mpi_core_changes.md)
- [test_application_changes](../../../../functions/tests/test_ci_utilities/TestCIUtilities/test_application_changes.md)