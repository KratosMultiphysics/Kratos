---
type: C++ Method
title: MapLocalToGlobal2D
resource: utilities/beam_math_utilities.hpp#L429-L441
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/quaternion/ToRotationMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void MapLocalToGlobal2D(const QuaternionType& rLocalToGlobalQuaternion, MatrixType& rMatrix)`

# Calls

- [ToRotationMatrix](../../../../functions/utilities/quaternion/ToRotationMatrix.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)