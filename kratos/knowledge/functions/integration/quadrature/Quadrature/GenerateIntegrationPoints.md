---
type: C++ Method
title: GenerateIntegrationPoints
resource: integration/quadrature.h#L125-L130
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/hexahedra_3d_20/Hexahedra3D20/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/hexahedra_3d_27/Hexahedra3D27/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/hexahedra_3d_8/Hexahedra3D8/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/hexahedra_interface_3d_8/HexahedraInterface3D8/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_2/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_3/Line2D3/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_4/Line2D4/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_5/Line2D5/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_3d_2/Line3D2/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_3d_3/Line3D3/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/prism_3d_15/Prism3D15/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/prism_3d_6/Prism3D6/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/prism_interface_3d_6/PrismInterface3D6/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/pyramid_3d_13/Pyramid3D13/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/pyramid_3d_5/Pyramid3D5/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_2d_4/Quadrilateral2D4/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_2d_8/Quadrilateral2D8/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_2d_9/Quadrilateral2D9/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_3d_4/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_3d_8/Quadrilateral3D8/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_3d_9/Quadrilateral3D9/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_interface_2d_4/QuadrilateralInterface2D4/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_interface_3d_4/QuadrilateralInterface3D4/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/sphere_3d_1/Sphere3D1/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/tetrahedra_3d_10/Tetrahedra3D10/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_10/Triangle2D10/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_15/Triangle2D15/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_3/Triangle2D3/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_6/Triangle2D6/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_6/Triangle3D6/AllIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/InsertGaussPointsExactSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/InsertGaussPointsApproxSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/GetIntegrationTriangle
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static const IntegrationPointsArrayType GenerateIntegrationPoints()`

# Called by

- [AllIntegrationPoints](../../../../functions/geometries/hexahedra_3d_20/Hexahedra3D20/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/hexahedra_3d_27/Hexahedra3D27/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/hexahedra_3d_8/Hexahedra3D8/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/hexahedra_interface_3d_8/HexahedraInterface3D8/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/line_2d_2/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/line_2d_3/Line2D3/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/line_2d_4/Line2D4/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/line_2d_5/Line2D5/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/line_3d_2/Line3D2/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/line_3d_3/Line3D3/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/prism_3d_15/Prism3D15/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/prism_3d_6/Prism3D6/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/prism_interface_3d_6/PrismInterface3D6/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/pyramid_3d_13/Pyramid3D13/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/pyramid_3d_5/Pyramid3D5/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/quadrilateral_2d_4/Quadrilateral2D4/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/quadrilateral_2d_8/Quadrilateral2D8/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/quadrilateral_2d_9/Quadrilateral2D9/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/quadrilateral_3d_4/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/quadrilateral_3d_8/Quadrilateral3D8/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/quadrilateral_3d_9/Quadrilateral3D9/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/quadrilateral_interface_2d_4/QuadrilateralInterface2D4/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/quadrilateral_interface_3d_4/QuadrilateralInterface3D4/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/sphere_3d_1/Sphere3D1/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/tetrahedra_3d_10/Tetrahedra3D10/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/triangle_2d_10/Triangle2D10/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/triangle_2d_15/Triangle2D15/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/triangle_2d_3/Triangle2D3/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/triangle_2d_6/Triangle2D6/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/triangle_3d_3/AllIntegrationPoints.md)
- [AllIntegrationPoints](../../../../functions/geometries/triangle_3d_6/Triangle3D6/AllIntegrationPoints.md)
- [InsertGaussPointsExactSurface](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/InsertGaussPointsExactSurface.md)
- [InsertGaussPointsApproxSurface](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/InsertGaussPointsApproxSurface.md)
- [GetIntegrationTriangle](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/GetIntegrationTriangle.md)