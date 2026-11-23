---
type: C++ Method
title: EigenVectors
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L2617-L2766
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void EigenVectors(const Matrix& A, Matrix& vectors, Vector& lambda, double zero_tolerance =1e-9, int max_iterations = 10)`

# Called by

- [CalcSignedDistancesToMoreThanThreeIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes.md)