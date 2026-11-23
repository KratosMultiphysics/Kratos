---
type: C++ Function
title: InverseCayleyTransform
resource: utilities/beam_math_utilities.hpp#L199-L220
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/InvertMatrix3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/beam_math_utilities/SkewSymmetricTensorToVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void InverseCayleyTransform(const MatrixType& rCayleyTensor, TVector3& rVector)`

# Calls

- [InvertMatrix3](../../../functions/utilities/math_utils/InvertMatrix3.md)
- [prod](../../../functions/includes/amatrix_interface/prod.md)
- [SkewSymmetricTensorToVector](../../../functions/utilities/beam_math_utilities/SkewSymmetricTensorToVector.md)