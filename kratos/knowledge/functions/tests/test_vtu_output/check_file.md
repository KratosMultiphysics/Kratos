---
type: Python Function
title: check_file
resource: tests/test_vtu_output.py#L116-L139
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/Check
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def check_file(output_file_name: str, reference_file_name: str):`

# Calls

- [SetString](../../../functions/sources/kratos_parameters/Parameters/SetString.md)

# Called by

- [Check](../../../functions/tests/test_vtu_output/TestVtuOutputBase/Check.md)