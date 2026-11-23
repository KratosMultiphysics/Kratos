---
type: C++ Method
title: KnotsU
resource: geometries/brep_surface.h#L344-L349
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateQuadraturePointGeometriesSBM
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/KnotRefinementU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CalculateEstimatedKnotLengthness
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometriesSbm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/ProjectionPointGlobalToLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_geometries_to_python/AddGeometriesToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateOuterInclinedRefinedTrimmedBrepSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateInnerRectRefinedTrimmedBrepSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_surface/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_volume_knot_refinement/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_geometries/TestGeometry/test_nurbs_surface_3d
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<double> KnotsU() const`

# Called by

- [CreateQuadraturePointGeometries](../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateQuadraturePointGeometries.md)
- [CreateQuadraturePointGeometriesSBM](../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateQuadraturePointGeometriesSBM.md)
- [KnotRefinementU](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/KnotRefinementU.md)
- [DegreeElevationU](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationU.md)
- [KnotRefinementU](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementU.md)
- [CalculateEstimatedKnotLengthness](../../../../functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CalculateEstimatedKnotLengthness.md)
- [CreateQuadraturePointGeometries](../../../../functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometries.md)
- [CreateQuadraturePointGeometriesSbm](../../../../functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometriesSbm.md)
- [ProjectionPointGlobalToLocalSpace](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/ProjectionPointGlobalToLocalSpace.md)
- [CreateQuadraturePointGeometries](../../../../functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/CreateQuadraturePointGeometries.md)
- [GetBrepSurfaceParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters.md)
- [AddGeometriesToPython](../../../../functions/python/add_geometries_to_python/AddGeometriesToPython.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_brep_surface/KRATOS_TEST_CASE_IN_SUITE.md)
- [GenerateOuterInclinedRefinedTrimmedBrepSurface](../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateOuterInclinedRefinedTrimmedBrepSurface.md)
- [GenerateInnerRectRefinedTrimmedBrepSurface](../../../../functions/tests/cpp_tests/geometries/test_brep_surface_trimmed/GenerateInnerRectRefinedTrimmedBrepSurface.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_surface/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_volume_knot_refinement/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_nurbs_surface_3d](../../../../functions/tests/test_geometries/TestGeometry/test_nurbs_surface_3d.md)