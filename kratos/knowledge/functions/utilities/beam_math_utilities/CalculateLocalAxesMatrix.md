---
type: C++ Function
title: CalculateLocalAxesMatrix
resource: utilities/beam_math_utilities.hpp#L686-L711
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/beam_math_utilities/CalculateLocalAxesVectors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/beam_math_utilities/CalculateLocalAxisVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void CalculateLocalAxesMatrix(const TVector3& rLocalX, MatrixType& rRotationMatrix)`

# Calls

- [CalculateLocalAxesVectors](../../../functions/utilities/beam_math_utilities/CalculateLocalAxesVectors.md)
- [CalculateLocalAxisVector](../../../functions/utilities/beam_math_utilities/CalculateLocalAxisVector.md)