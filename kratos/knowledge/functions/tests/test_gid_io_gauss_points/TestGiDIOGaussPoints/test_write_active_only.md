---
type: Python Method
title: test_write_active_only
resource: tests/test_gid_io_gauss_points.py#L119-L132
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/deactivateSome
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/initializeOutputFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/writeResults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/finalizeOutputFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/outputMatchesReferenceSolution
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_write_active_only(self):`

# Calls

- [deactivateSome](../../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/deactivateSome.md)
- [initializeOutputFile](../../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/initializeOutputFile.md)
- [writeResults](../../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/writeResults.md)
- [finalizeOutputFile](../../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/finalizeOutputFile.md)
- [outputMatchesReferenceSolution](../../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/outputMatchesReferenceSolution.md)