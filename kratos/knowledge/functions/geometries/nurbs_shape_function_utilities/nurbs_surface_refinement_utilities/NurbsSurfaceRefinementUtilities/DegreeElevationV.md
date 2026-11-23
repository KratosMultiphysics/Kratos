---
type: C++ Method
title: DegreeElevationV
resource: geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities.cpp#L727-L1107
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
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

`void NurbsSurfaceRefinementUtilities::DegreeElevationV( NurbsSurfaceGeometryType& rGeometry, SizeType& rDegreeVToElevate, PointerVector<NodeType>& rPointsRefined, Vector& rKnotsVRefined, Vector& rWeightsRefined)`

# Calls

- [KnotsV](../../../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [fill](../../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [GetBinomCoefficient](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetBinomCoefficient.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/geometries/test_nurbs_surface/KRATOS_TEST_CASE_IN_SUITE.md)