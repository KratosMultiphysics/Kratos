---
type: C++ Function
title: GaussSeidelEigenSystem
resource: utilities/math_utils.h#L1682-L1790
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_metric_calculator/GeometryMetricCalculator<TDim,TNumNodes>/CalculateMetricTensor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/EigenSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/MatrixSquareRoot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/plane_approximation_utility/PlaneApproximationUtility/ComputePlaneNormal
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline bool GaussSeidelEigenSystem( const TMatrixType1& rA, TMatrixType2& rEigenVectorsMatrix, TMatrixType2& rEigenValuesMatrix, const double Tolerance = 1.0e-18, const SizeType MaxIterations = 20 )`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [CalculateMetricTensor](../../../functions/utilities/geometry_metric_calculator/GeometryMetricCalculator<TDim,TNumNodes>/CalculateMetricTensor.md)
- [EigenSystem](../../../functions/utilities/math_utils/EigenSystem.md)
- [MatrixSquareRoot](../../../functions/utilities/math_utils/MatrixSquareRoot.md)
- [ComputePlaneNormal](../../../functions/utilities/plane_approximation_utility/PlaneApproximationUtility/ComputePlaneNormal.md)