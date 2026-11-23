---
type: C++ Function
title: EigenSystem
resource: utilities/math_utils.h#L1806-L1824
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/GaussSeidelEigenSystem
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline bool EigenSystem( const BoundedMatrix<double, TDim, TDim>& rA, BoundedMatrix<double, TDim, TDim>& rEigenVectorsMatrix, BoundedMatrix<double, TDim, TDim>& rEigenValuesMatrix, const double Tolerance = 1.0e-18, const SizeType MaxIterations = 20 )`

# Calls

- [GaussSeidelEigenSystem](../../../functions/utilities/math_utils/GaussSeidelEigenSystem.md)