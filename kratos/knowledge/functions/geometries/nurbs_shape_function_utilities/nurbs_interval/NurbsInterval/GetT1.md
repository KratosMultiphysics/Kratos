---
type: C++ Method
title: GetT1
resource: geometries/nurbs_shape_function_utilities/nurbs_interval.h#L70-L73
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/brep_curve_on_surface/BrepCurveOnSurface/SpansLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface_on_volume/BrepSurfaceOnVolume/SpansLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/SpansLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/SpansLocalSpaceSBM
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_curve_axis_intersection/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double GetT1() const`

# Called by

- [SpansLocalSpace](../../../../../functions/geometries/brep_curve_on_surface/BrepCurveOnSurface/SpansLocalSpace.md)
- [SpansLocalSpace](../../../../../functions/geometries/brep_surface_on_volume/BrepSurfaceOnVolume/SpansLocalSpace.md)
- [SpansLocalSpace](../../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/SpansLocalSpace.md)
- [SpansLocalSpaceSBM](../../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/SpansLocalSpaceSBM.md)
- [CreateIntegrationPoints](../../../../../functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CreateIntegrationPoints.md)
- [CreateIntegrationPoints](../../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CreateIntegrationPoints.md)
- [GetBoundaryLoopParameters](../../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/utilities/test_curve_axis_intersection/KRATOS_TEST_CASE_IN_SUITE.md)