---
type: C++ Function
title: Norm3
resource: utilities/math_utils.h#L786-L791
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/hexahedra_3d_8/Hexahedra3D8/AverageEdgeLength
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/hexahedra_interface_3d_8/HexahedraInterface3D8/HexahedraInterface3D8
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/hexahedra_interface_3d_8/HexahedraInterface3D8/DeterminantOfJacobian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/prism_interface_3d_6/PrismInterface3D6/Area
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_interface_3d_4/QuadrilateralInterface3D4/QuadrilateralInterface3D4
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_interface_3d_4/QuadrilateralInterface3D4/Length
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/tetrahedra_3d_10/Tetrahedra3D10/FacesArePlanar
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/AverageEdgeLength
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_3/Triangle2D3/Semiperimeter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_3/Triangle2D3/AverageEdgeLength
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_3/Triangle2D3/Circumradius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_3/Triangle2D3/Inradius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_2d_3/Triangle2D3/InradiusToCircumradiusQuality
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/Area
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/AverageEdgeLength
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/Circumradius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/Inradius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/InradiusToCircumradiusQuality
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_3/AreaToEdgeLengthSquareRatio
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/triangle_3d_6/Triangle3D6/EdgesAreStraight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeSquaredArea
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_quadrilateral_3d_9/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/intersection_utilities/IntersectionUtilities/ComputeTriangleLineIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/HeronCheck
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double Norm3(const TVectorType& rA)`

# Called by

- [AverageEdgeLength](../../../functions/geometries/hexahedra_3d_8/Hexahedra3D8/AverageEdgeLength.md)
- [HexahedraInterface3D8](../../../functions/geometries/hexahedra_interface_3d_8/HexahedraInterface3D8/HexahedraInterface3D8.md)
- [DeterminantOfJacobian](../../../functions/geometries/hexahedra_interface_3d_8/HexahedraInterface3D8/DeterminantOfJacobian.md)
- [Area](../../../functions/geometries/prism_interface_3d_6/PrismInterface3D6/Area.md)
- [QuadrilateralInterface3D4](../../../functions/geometries/quadrilateral_interface_3d_4/QuadrilateralInterface3D4/QuadrilateralInterface3D4.md)
- [Length](../../../functions/geometries/quadrilateral_interface_3d_4/QuadrilateralInterface3D4/Length.md)
- [FacesArePlanar](../../../functions/geometries/tetrahedra_3d_10/Tetrahedra3D10/FacesArePlanar.md)
- [AverageEdgeLength](../../../functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/AverageEdgeLength.md)
- [Semiperimeter](../../../functions/geometries/triangle_2d_3/Triangle2D3/Semiperimeter.md)
- [AverageEdgeLength](../../../functions/geometries/triangle_2d_3/Triangle2D3/AverageEdgeLength.md)
- [Circumradius](../../../functions/geometries/triangle_2d_3/Triangle2D3/Circumradius.md)
- [Inradius](../../../functions/geometries/triangle_2d_3/Triangle2D3/Inradius.md)
- [InradiusToCircumradiusQuality](../../../functions/geometries/triangle_2d_3/Triangle2D3/InradiusToCircumradiusQuality.md)
- [Area](../../../functions/geometries/triangle_3d_3/Area.md)
- [AverageEdgeLength](../../../functions/geometries/triangle_3d_3/AverageEdgeLength.md)
- [Circumradius](../../../functions/geometries/triangle_3d_3/Circumradius.md)
- [Inradius](../../../functions/geometries/triangle_3d_3/Inradius.md)
- [InradiusToCircumradiusQuality](../../../functions/geometries/triangle_3d_3/InradiusToCircumradiusQuality.md)
- [AreaToEdgeLengthSquareRatio](../../../functions/geometries/triangle_3d_3/AreaToEdgeLengthSquareRatio.md)
- [EdgesAreStraight](../../../functions/geometries/triangle_3d_6/Triangle3D6/EdgesAreStraight.md)
- [ComputeSquaredArea](../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeSquaredArea.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_quadrilateral_3d_9/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [ComputeTriangleLineIntersection](../../../functions/utilities/intersection_utilities/IntersectionUtilities/ComputeTriangleLineIntersection.md)
- [HeronCheck](../../../functions/utilities/mortar_utilities/HeronCheck.md)