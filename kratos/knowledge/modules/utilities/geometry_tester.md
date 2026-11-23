---
type: C++ Module
title: geometry_tester
resource: utilities/geometry_tester.cpp#L1-L1144
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/spaces-ublas-space-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-geometry-tester-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-element-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-geometry-data-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-geometry-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-geometry-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-2d-3-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-2d-6-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-2d-4-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-2d-9-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-interface-2d-4-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-tetrahedra-3d-4-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-tetrahedra-3d-10-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-hexahedra-3d-8-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-hexahedra-3d-20-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-hexahedra-3d-27-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-hexahedra-interface-3d-8-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-prism-3d-6-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-prism-3d-15-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-prism-interface-3d-6-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [RunTest](../../functions/utilities/geometry_tester/GeometryTesterUtility/RunTest.md)
- [StreamTestTetrahedra3D4N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTetrahedra3D4N.md)
- [StreamTestTetrahedra3D10N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTetrahedra3D10N.md)
- [StreamTestTriangle2D3N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTriangle2D3N.md)
- [StreamTestTriangle2D6N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTriangle2D6N.md)
- [StreamTestQuadrilateral2D4N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateral2D4N.md)
- [StreamTestQuadrilateral2D9N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateral2D9N.md)
- [StreamTestQuadrilateralInterface2D4N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateralInterface2D4N.md)
- [StreamTestHexahedra3D8N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D8N.md)
- [StreamTestHexahedra3D20N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D20N.md)
- [StreamTestHexahedra3D27N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D27N.md)
- [StreamTestHexahedraInterface3D8N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedraInterface3D8N.md)
- [StreamTestPrism3D6N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestPrism3D6N.md)
- [StreamTestPrism3D15N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestPrism3D15N.md)
- [StreamTestPrismInterface3D6N](../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestPrismInterface3D6N.md)
- [GenerateNodes](../../functions/utilities/geometry_tester/GeometryTesterUtility/GenerateNodes.md)
- [VerifyAreaByIntegration](../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyAreaByIntegration.md)
- [VerifyStrainExactness](../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyStrainExactness.md)
- [CalculateB](../../functions/utilities/geometry_tester/GeometryTesterUtility/CalculateB.md)
- [VerifyIsInside](../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyIsInside.md)
- [VerfiyShapeFunctionsValues](../../functions/utilities/geometry_tester/GeometryTesterUtility/VerfiyShapeFunctionsValues.md)
- [VerifyShapeFunctionsSecondDerivativesValues](../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesValues.md)
- [VerifyShapeFunctionsSecondDerivativesInterpolation](../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesInterpolation.md)
- [GetIntegrationName](../../functions/utilities/geometry_tester/GeometryTesterUtility/GetIntegrationName.md)
- [GetGeometryName](../../functions/utilities/geometry_tester/GeometryTesterUtility/GetGeometryName.md)

# Imports

- `"spaces/ublas_space.h"`
- `"utilities/geometry_tester.h"`
- `"includes/element.h"`
- `"geometries/geometry_data.h"`
- `"utilities/geometry_utilities.h"`
- `"geometries/geometry.h"`
- `"geometries/triangle_2d_3.h"`
- `"geometries/triangle_2d_6.h"`
- `"geometries/quadrilateral_2d_4.h"`
- `"geometries/quadrilateral_2d_9.h"`
- `"geometries/quadrilateral_interface_2d_4.h"`
- `"geometries/tetrahedra_3d_4.h"`
- `"geometries/tetrahedra_3d_10.h"`
- `"geometries/hexahedra_3d_8.h"`
- `"geometries/hexahedra_3d_20.h"`
- `"geometries/hexahedra_3d_27.h"`
- `"geometries/hexahedra_interface_3d_8.h"`
- `"geometries/prism_3d_6.h"`
- `"geometries/prism_3d_15.h"`
- `"geometries/prism_interface_3d_6.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)