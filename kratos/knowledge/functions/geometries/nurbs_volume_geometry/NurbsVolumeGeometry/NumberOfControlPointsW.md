---
type: C++ Method
title: NumberOfControlPointsW
resource: geometries/nurbs_volume_geometry.h#L419-L422
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfKnotsW
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometriesSbm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/PointsNumberInDirection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/GlobalSpaceDerivatives
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_volume_knot_refinement/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SizeType NumberOfControlPointsW() const`

# Calls

- [NumberOfKnotsW](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfKnotsW.md)

# Called by

- [KnotRefinementU](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementU.md)
- [KnotRefinementV](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementV.md)
- [KnotRefinementW](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementW.md)
- [CreateQuadraturePointGeometries](../../../../functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometries.md)
- [CreateQuadraturePointGeometriesSbm](../../../../functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometriesSbm.md)
- [PointsNumberInDirection](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/PointsNumberInDirection.md)
- [CreateQuadraturePointGeometries](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CreateQuadraturePointGeometries.md)
- [GlobalSpaceDerivatives](../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/GlobalSpaceDerivatives.md)
- [CreateQuadraturePointGeometries](../../../../functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/CreateQuadraturePointGeometries.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_volume_knot_refinement/KRATOS_TEST_CASE_IN_SUITE.md)