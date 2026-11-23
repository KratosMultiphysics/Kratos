---
type: C++ Function
title: ExplicitSolveUpdate
resource: processes/flux_corrected_transport_convection_process.h#L481-L492
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/flux_corrected_transport_convection_process/ExplicitSolve
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ExplicitSolveUpdate( std::vector<double>& rUpdate, const std::vector<std::vector<double>>& rResiduals, const double DeltaTime, const TSolveCoeffsIterType& rSolveCoeffBegin, const TSolveCoeffsIterType& rSolveCoeffEnd, const IndexType MaxIteration)`

# Calls

- [ExplicitSolve](../../../functions/processes/flux_corrected_transport_convection_process/ExplicitSolve.md)

# Called by

- [SolveSubStep](../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep.md)