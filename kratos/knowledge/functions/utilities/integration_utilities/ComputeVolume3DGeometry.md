---
type: C++ Function
title: ComputeVolume3DGeometry
resource: utilities/integration_utilities.h#L168-L180
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/Det3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/hexahedra_3d_20/Hexahedra3D20/Volume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/hexahedra_3d_27/Hexahedra3D27/Volume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/hexahedra_3d_8/Hexahedra3D8/Volume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/prism_3d_15/Prism3D15/Volume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/prism_3d_6/Prism3D6/Volume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/pyramid_3d_13/Pyramid3D13/Volume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/pyramid_3d_5/Pyramid3D5/Volume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/tetrahedra_3d_10/Tetrahedra3D10/Volume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_integration_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double ComputeVolume3DGeometry(const Geometry<TPointType>& rGeometry)`

# Calls

- [GetDefaultIntegrationMethod](../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [Det3](../../../functions/utilities/math_utils/Det3.md)
- [Weight](../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [Volume](../../../functions/geometries/hexahedra_3d_20/Hexahedra3D20/Volume.md)
- [Volume](../../../functions/geometries/hexahedra_3d_27/Hexahedra3D27/Volume.md)
- [Volume](../../../functions/geometries/hexahedra_3d_8/Hexahedra3D8/Volume.md)
- [Volume](../../../functions/geometries/prism_3d_15/Prism3D15/Volume.md)
- [Volume](../../../functions/geometries/prism_3d_6/Prism3D6/Volume.md)
- [Volume](../../../functions/geometries/pyramid_3d_13/Pyramid3D13/Volume.md)
- [Volume](../../../functions/geometries/pyramid_3d_5/Pyramid3D5/Volume.md)
- [Volume](../../../functions/geometries/tetrahedra_3d_10/Tetrahedra3D10/Volume.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_integration_utilities/KRATOS_TEST_CASE_IN_SUITE.md)