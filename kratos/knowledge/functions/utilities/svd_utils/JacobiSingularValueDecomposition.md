---
type: C++ Function
title: JacobiSingularValueDecomposition
resource: utilities/svd_utils.h#L178-L248
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/svd_utils/ComputeJacobiIterations
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_svd_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/qef_utility/QuadraticErrorFunction/ComputeQuadraticErrorFunctionPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/svd_utils/SingularValueDecomposition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline IndexType JacobiSingularValueDecomposition( const TInputMatrix& rInputMatrix, TUmatrix& rUMatrix, TSMatrix& rSMatrix, TVMatrix& rVMatrix, const TDataType Tolerance = ZeroTolerance, const IndexType MaxIter = 200 )`

# Calls

- [ComputeJacobiIterations](../../../functions/utilities/svd_utils/ComputeJacobiIterations.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_svd_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [ComputeQuadraticErrorFunctionPoint](../../../functions/utilities/qef_utility/QuadraticErrorFunction/ComputeQuadraticErrorFunctionPoint.md)
- [SingularValueDecomposition](../../../functions/utilities/svd_utils/SingularValueDecomposition.md)