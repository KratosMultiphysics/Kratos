---
type: C++ Method
title: EvaluateLimiter
resource: processes/flux_corrected_transport_convection_process.h#L802-L919
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/GetLumpedMassMatrixDiagonal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeData/GetAntidiffusiveEdgeContribution
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EvaluateLimiter(const double DeltaTime)`

# Calls

- [GetLumpedMassMatrixDiagonal](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/GetLumpedMassMatrixDiagonal.md)
- [GetAntidiffusiveEdgeContribution](../../../../functions/containers/edge_based_data_structure/EdgeData/GetAntidiffusiveEdgeContribution.md)