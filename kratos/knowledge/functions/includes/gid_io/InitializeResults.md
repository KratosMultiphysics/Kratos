---
type: C++ Function
title: InitializeResults
resource: includes/gid_io.h#L566-L613
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/gid_gauss_point_container/GidGaussPointsContainer/WriteGaussPoints
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output/GiDOutput/_initialize_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__initialize_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/initializeOutputFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void InitializeResults( const double name, const MeshType& rThisMesh )`

# Calls

- [WriteGaussPoints](../../../functions/includes/gid_gauss_point_container/GidGaussPointsContainer/WriteGaussPoints.md)

# Called by

- [_initialize_results](../../../functions/python_scripts/gid_output/GiDOutput/_initialize_results.md)
- [__initialize_results](../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__initialize_results.md)
- [initializeOutputFile](../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/initializeOutputFile.md)