---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_nurbs_surface.cpp#L304-L332
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_nurbs_surface/GenerateReferencePointSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_surface/GenerateReferencePieceOfCylinderNurbsSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_surface/GenerateReferenceNodeSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_surface/GenerateReferenceQuarterSphereGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_surface/GenerateReferenceArbitraryNodeSurface
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(NurbsSurfacePoint, KratosCoreNurbsGeometriesFastSuite)`

# Calls

- [GenerateReferencePointSurface](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_surface/GenerateReferencePointSurface.md)
- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)
- [GenerateReferencePieceOfCylinderNurbsSurface](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_surface/GenerateReferencePieceOfCylinderNurbsSurface.md)
- [norm_2](../../../../../functions/includes/amatrix_interface/norm_2.md)
- [GenerateReferenceNodeSurface](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_surface/GenerateReferenceNodeSurface.md)
- [Weight](../../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [KnotsV](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [KnotsU](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [DegreeElevationU](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationU.md)
- [DegreeElevationV](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationV.md)
- [GenerateReferenceQuarterSphereGeometry](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_surface/GenerateReferenceQuarterSphereGeometry.md)
- [GenerateReferenceArbitraryNodeSurface](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_surface/GenerateReferenceArbitraryNodeSurface.md)