---
type: Python Module
title: test_gid_io_gauss_points
resource: tests/test_gid_io_gauss_points.py#L1-L153
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratosunittest-as-unittest
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratos-utilities-as-kratos-utils
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-fluiddynamicsapplication
    resolved_by: tree-sitter
    confidence: exact
  - target: external/filecmp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/os
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [WorkFolderScope](../../classes/tests/test_gid_io_gauss_points/WorkFolderScope.md)
- [__init__](../../functions/tests/test_gid_io_gauss_points/WorkFolderScope/__init__.md)
- [__enter__](../../functions/tests/test_gid_io_gauss_points/WorkFolderScope/__enter__.md)
- [__exit__](../../functions/tests/test_gid_io_gauss_points/WorkFolderScope/__exit__.md)
- [TestGiDIOGaussPoints](../../classes/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints.md)
- [setUp](../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/setUp.md)
- [tearDown](../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/tearDown.md)
- [setModelPart](../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/setModelPart.md)
- [deactivateSome](../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/deactivateSome.md)
- [initializeOutputFile](../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/initializeOutputFile.md)
- [writeResults](../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/writeResults.md)
- [finalizeOutputFile](../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/finalizeOutputFile.md)
- [outputMatchesReferenceSolution](../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/outputMatchesReferenceSolution.md)
- [test_write_active_only](../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/test_write_active_only.md)
- [test_write_dynamic_deactivation](../../functions/tests/test_gid_io_gauss_points/TestGiDIOGaussPoints/test_write_dynamic_deactivation.md)

# Imports

- `KratosMultiphysics`
- `KratosMultiphysics.KratosUnittest as UnitTest`
- `KratosMultiphysics.kratos_utilities as kratos_utils`
- `KratosMultiphysics.FluidDynamicsApplication`
- `filecmp`
- `os`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)