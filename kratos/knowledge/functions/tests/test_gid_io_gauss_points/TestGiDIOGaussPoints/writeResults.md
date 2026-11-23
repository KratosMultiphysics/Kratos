---
type: Python Method
title: writeResults
resource: tests/test_gid_io_gauss_points.py#L103-L107
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/gid_io/PrintFlagsOnGaussPoints
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

`def writeResults(self,label):`

# Calls

- [PrintFlagsOnGaussPoints](../../../../functions/includes/gid_io/PrintFlagsOnGaussPoints.md)

# Called by

- [test_write_active_only](../../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/test_write_active_only.md)
- [test_write_dynamic_deactivation](../../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/test_write_dynamic_deactivation.md)