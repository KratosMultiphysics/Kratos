---
type: C++ Function
title: Det3
resource: utilities/math_utils.h#L578-L588
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/geometries/test_hexahedra_3d_27/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/integration_utilities/ComputeVolume3DGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double Det3(const TMatrixType& rA)`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_hexahedra_3d_27/KRATOS_TEST_CASE_IN_SUITE.md)
- [ComputeVolume3DGeometry](../../../functions/utilities/integration_utilities/ComputeVolume3DGeometry.md)