---
type: C++ Method
title: ComputeTessellation
resource: utilities/tessellation_utilities/curve_tessellation.h#L162-L262
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputeGlobalCoordinates
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector/emplace_back
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_set/PointerVectorSet/pop_back
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/DistanceToLine
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeBoundaryTessellation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_curve_tessellation/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/Tessellate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static TessellationType ComputeTessellation( const GeometryType& rGeometry, const int PolynomialDegree, const NurbsInterval DomainInterval, const std::vector<NurbsInterval>& rKnotSpanIntervals, const double Tolerance, bool ToSurfaceParameter = false )`

# Calls

- [GetT0](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/GetT0.md)
- [ComputeGlobalCoordinates](../../../../../functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputeGlobalCoordinates.md)
- [emplace_back](../../../../../functions/containers/pointer_vector/emplace_back.md)
- [pop_back](../../../../../functions/containers/pointer_vector_set/PointerVectorSet/pop_back.md)
- [DistanceToLine](../../../../../functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/DistanceToLine.md)

# Called by

- [ComputeBoundaryTessellation](../../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeBoundaryTessellation.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/utilities/test_curve_tessellation/KRATOS_TEST_CASE_IN_SUITE.md)
- [Tessellate](../../../../../functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/Tessellate.md)