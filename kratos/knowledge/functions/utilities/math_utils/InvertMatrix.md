---
type: C++ Function
title: InvertMatrix
resource: utilities/math_utils.h#L213-L243
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/InvertMatrix2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/InvertMatrix3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/InvertMatrix4
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/CheckConditionNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline BoundedMatrix<double, TDim, TDim> InvertMatrix( const BoundedMatrix<double, TDim, TDim>& rInputMatrix, double& rInputMatrixDet, const double Tolerance = ZeroTolerance )`

# Calls

- [InvertMatrix2](../../../functions/utilities/math_utils/InvertMatrix2.md)
- [InvertMatrix3](../../../functions/utilities/math_utils/InvertMatrix3.md)
- [InvertMatrix4](../../../functions/utilities/math_utils/InvertMatrix4.md)
- [CheckConditionNumber](../../../functions/utilities/math_utils/CheckConditionNumber.md)