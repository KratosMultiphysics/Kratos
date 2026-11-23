---
type: Python Method
title: finalizeOutputFile
resource: tests/test_gid_io_gauss_points.py#L110-L111
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/gid_io/FinalizeResults
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/test_write_active_only
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/test_write_dynamic_deactivation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def finalizeOutputFile(self):`

# Calls

- [FinalizeResults](../../../../functions/includes/gid_io/FinalizeResults.md)

# Called by

- [test_write_active_only](../../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/test_write_active_only.md)
- [test_write_dynamic_deactivation](../../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/test_write_dynamic_deactivation.md)