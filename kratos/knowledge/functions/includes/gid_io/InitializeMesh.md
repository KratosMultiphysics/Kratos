---
type: C++ Function
title: InitializeMesh
resource: includes/gid_io.h#L1054-L1101
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/gid_output/GiDOutput/_write_mesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__write_mesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/unv_output_process/UnvOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_tetrahedra_edge_swapping_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/initializeOutputFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InitializeMesh( const double name )`

# Called by

- [_write_mesh](../../../functions/python_scripts/gid_output/GiDOutput/_write_mesh.md)
- [__write_mesh](../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_mesh.md)
- [__init__](../../../functions/python_scripts/unv_output_process/UnvOutputProcess/__init__.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/processes/test_tetrahedra_edge_swapping_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [initializeOutputFile](../../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/initializeOutputFile.md)
- [TestIODebug](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug.md)