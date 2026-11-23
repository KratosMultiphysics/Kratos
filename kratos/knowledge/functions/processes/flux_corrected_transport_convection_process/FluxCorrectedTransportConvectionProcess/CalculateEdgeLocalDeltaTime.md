---
type: C++ Method
title: CalculateEdgeLocalDeltaTime
resource: processes/flux_corrected_transport_convection_process.h#L359-L368
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateSubStepDeltaTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateResidual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateAntidiffusiveEdgeContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double CalculateEdgeLocalDeltaTime( const double Length, const double NormVelI, const double NormVelJ)`

# Called by

- [CalculateSubStepDeltaTime](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateSubStepDeltaTime.md)
- [CalculateResidual](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateResidual.md)
- [CalculateAntidiffusiveEdgeContributions](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateAntidiffusiveEdgeContributions.md)