---
type: Python Function
title: changed_files
resource: python_scripts/testing/ci_utilities.py#L26-L28
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/testing/ci_utilities/check_valid_environment_configuration_exists
    resolved_by: tree-sitter
    confidence: exact
  called_by:
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
  - target: functions/python_scripts/testing/ci_utilities/print_ci_information
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def changed_files() -> List[Path]:`

# Calls

- [check_valid_environment_configuration_exists](../../../../functions/python_scripts/testing/ci_utilities/check_valid_environment_configuration_exists.md)

# Called by

- [get_changed_applications](../../../../functions/python_scripts/testing/ci_utilities/get_changed_applications.md)
- [is_core_changed](../../../../functions/python_scripts/testing/ci_utilities/is_core_changed.md)
- [is_mpi_core_changed](../../../../functions/python_scripts/testing/ci_utilities/is_mpi_core_changed.md)
- [get_changed_files_extensions](../../../../functions/python_scripts/testing/ci_utilities/get_changed_files_extensions.md)
- [print_ci_information](../../../../functions/python_scripts/testing/ci_utilities/print_ci_information.md)