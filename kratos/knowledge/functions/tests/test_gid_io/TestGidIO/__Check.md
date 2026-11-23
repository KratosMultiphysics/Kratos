---
type: Python Method
title: __Check
resource: tests/test_gid_io.py#L61-L82
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_gid_io/TestGidIO/test_gid_io_all
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io/TestGidIO/test_gid_io_deactivation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io/TestGidIO/test_gid_io_results
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __Check(self,output_file,reference_file): ## Settings string in json format`

# Calls

- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)

# Called by

- [test_gid_io_all](../../../../functions/tests/test_gid_io/TestGidIO/test_gid_io_all.md)
- [test_gid_io_deactivation](../../../../functions/tests/test_gid_io/TestGidIO/test_gid_io_deactivation.md)
- [test_gid_io_results](../../../../functions/tests/test_gid_io/TestGidIO/test_gid_io_results.md)