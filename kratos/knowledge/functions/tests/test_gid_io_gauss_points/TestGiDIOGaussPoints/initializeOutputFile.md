---
type: Python Method
title: initializeOutputFile
resource: tests/test_gid_io_gauss_points.py#L89-L101
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/gid_io/GidIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/InitializeMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/FinalizeMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/InitializeResults
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

`def initializeOutputFile(self):`

# Calls

- [GidIO](../../../../functions/includes/gid_io/GidIO.md)
- [InitializeMesh](../../../../functions/includes/gid_io/InitializeMesh.md)
- [FinalizeMesh](../../../../functions/includes/gid_io/FinalizeMesh.md)
- [InitializeResults](../../../../functions/includes/gid_io/InitializeResults.md)

# Called by

- [test_write_active_only](../../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/test_write_active_only.md)
- [test_write_dynamic_deactivation](../../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/test_write_dynamic_deactivation.md)