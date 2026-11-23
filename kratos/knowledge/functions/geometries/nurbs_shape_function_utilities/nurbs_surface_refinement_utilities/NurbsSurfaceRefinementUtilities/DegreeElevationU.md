---
type: C++ Method
title: DegreeElevationU
resource: geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities.cpp#L345-L725
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetBinomCoefficient
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/geometries/test_nurbs_surface/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NurbsSurfaceRefinementUtilities::DegreeElevationU( NurbsSurfaceGeometryType& rGeometry, SizeType& rDegreeUToElevate, PointerVector<NodeType>& rPointsRefined, Vector& rKnotsURefined, Vector& rWeightsRefined)`

# Calls

- [KnotsU](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [fill](../../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [GetBinomCoefficient](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetBinomCoefficient.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_surface/KRATOS_TEST_CASE_IN_SUITE.md)