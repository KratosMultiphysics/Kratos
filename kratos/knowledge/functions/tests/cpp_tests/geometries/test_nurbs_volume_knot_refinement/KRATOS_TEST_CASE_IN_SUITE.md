---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/geometries/test_nurbs_volume_knot_refinement.cpp#L97-L143
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/geometries/test_nurbs_volume_knot_refinement/GenerateVolume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementW
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(NurbsVolumeKnotRefinementU, KratosCoreNurbsGeometriesFastSuite)`

# Calls

- [GenerateVolume](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_volume_knot_refinement/GenerateVolume.md)
- [KnotsV](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [KnotsU](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [NumberOfControlPointsW](../../../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/NumberOfControlPointsW.md)
- [KnotRefinementW](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_volume_refinement_utilities/NurbsVolumeRefinementUtilities/KnotRefinementW.md)