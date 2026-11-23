---
type: C++ Method
title: SolveSubStep
resource: processes/flux_corrected_transport_convection_process.h#L370-L449
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/GetButcherTableauVariant
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/SubstepCount
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetIntegrationTheta
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateResidual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetMatrixRow
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/ExplicitSolveSolution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateAntidiffusiveEdgeContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/ExplicitSolveUpdate
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SolveSubStep(const double DeltaTime)`

# Calls

- [GetButcherTableauVariant](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/GetButcherTableauVariant.md)
- [visit](../../../../functions/python/variable_indexing_python/visit.md)
- [SubstepCount](../../../../functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/SubstepCount.md)
- [GetIntegrationTheta](../../../../functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetIntegrationTheta.md)
- [CalculateResidual](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateResidual.md)
- [GetMatrixRow](../../../../functions/solving_strategies/strategies/butcher_tableau/ButcherTableau/GetMatrixRow.md)
- [ExplicitSolveSolution](../../../../functions/processes/flux_corrected_transport_convection_process/ExplicitSolveSolution.md)
- [CalculateAntidiffusiveEdgeContributions](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateAntidiffusiveEdgeContributions.md)
- [ExplicitSolveUpdate](../../../../functions/processes/flux_corrected_transport_convection_process/ExplicitSolveUpdate.md)

# Called by

- [Execute](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/Execute.md)