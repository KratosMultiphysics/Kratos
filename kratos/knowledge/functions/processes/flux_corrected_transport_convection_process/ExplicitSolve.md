---
type: C++ Function
title: ExplicitSolve
resource: processes/flux_corrected_transport_convection_process.h#L495-L608
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/GetLumpedMassMatrixDiagonal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeData/GetOffDiagonalConsistentMass
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/ExplicitSolveSolution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/ExplicitSolveUpdate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ExplicitSolve( std::vector<double>& rSolution, const std::vector<std::vector<double>>& rResiduals, const double DeltaTime, const TSolveCoeffsIterType& rSolveCoeffBegin, const TSolveCoeffsIterType& rSolveCoeffEnd, const std::function<void(std::vector<double>&, const IndexType, const double)>& rUpdater, const IndexType MaxIteration)`

# Calls

- [GetLumpedMassMatrixDiagonal](../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/GetLumpedMassMatrixDiagonal.md)
- [GetOffDiagonalConsistentMass](../../../functions/containers/edge_based_data_structure/EdgeData/GetOffDiagonalConsistentMass.md)

# Called by

- [ExplicitSolveSolution](../../../functions/processes/flux_corrected_transport_convection_process/ExplicitSolveSolution.md)
- [ExplicitSolveUpdate](../../../functions/processes/flux_corrected_transport_convection_process/ExplicitSolveUpdate.md)