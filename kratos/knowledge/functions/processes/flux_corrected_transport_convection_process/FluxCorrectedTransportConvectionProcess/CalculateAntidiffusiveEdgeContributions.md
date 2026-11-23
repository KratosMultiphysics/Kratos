---
type: C++ Method
title: CalculateAntidiffusiveEdgeContributions
resource: processes/flux_corrected_transport_convection_process.h#L749-L800
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateEdgeLocalDeltaTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeData/GetOffDiagonalConsistentMass
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeData/SetAntidiffusiveEdgeContribution
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateAntidiffusiveEdgeContributions(const double DeltaTime)`

# Calls

- [CalculateEdgeLocalDeltaTime](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateEdgeLocalDeltaTime.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [GetOffDiagonalConsistentMass](../../../../functions/containers/edge_based_data_structure/EdgeData/GetOffDiagonalConsistentMass.md)
- [SetAntidiffusiveEdgeContribution](../../../../functions/containers/edge_based_data_structure/EdgeData/SetAntidiffusiveEdgeContribution.md)

# Called by

- [SolveSubStep](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep.md)