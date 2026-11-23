---
type: Python Function
title: ci_applications
resource: python_scripts/testing/ci_utilities.py#L31-L40
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/testing/ci_utilities/check_valid_environment_configuration_exists
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/testing/ci_utilities/print_ci_information
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ci_applications() -> List[str]:`

# Calls

- [check_valid_environment_configuration_exists](../../../../functions/python_scripts/testing/ci_utilities/check_valid_environment_configuration_exists.md)

# Called by

- [print_ci_information](../../../../functions/python_scripts/testing/ci_utilities/print_ci_information.md)