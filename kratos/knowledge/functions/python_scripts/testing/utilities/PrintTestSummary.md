---
type: Python Function
title: PrintTestSummary
resource: python_scripts/testing/utilities.py#L96-L104
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/testing/run_cpp_mpi_tests/main
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/run_cpp_tests/main
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/run_python_mpi_tests/main
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/run_python_tests/main
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def PrintTestSummary(exit_codes):`

# Called by

- [main](../../../../functions/python_scripts/testing/run_cpp_mpi_tests/main.md)
- [main](../../../../functions/python_scripts/testing/run_cpp_tests/main.md)
- [main](../../../../functions/python_scripts/testing/run_python_mpi_tests/main.md)
- [main](../../../../functions/python_scripts/testing/run_python_tests/main.md)