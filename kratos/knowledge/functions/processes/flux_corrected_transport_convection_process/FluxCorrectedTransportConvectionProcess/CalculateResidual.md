---
type: C++ Method
title: CalculateResidual
resource: processes/flux_corrected_transport_convection_process.h#L610-L747
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/edge_based_data_structure/EdgeData/IsBoundary
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateEdgeLocalDeltaTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeData/GetOffDiagonalConsistentMass
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateResidual( std::vector<std::vector<double>>& rResiduals, const std::vector<double>& rSolutionVector, const IndexType ResidualColumn, const double DiffusionConstant, const double DeltaTime)`

# Calls

- [IsBoundary](../../../../functions/containers/edge_based_data_structure/EdgeData/IsBoundary.md)
- [CalculateEdgeLocalDeltaTime](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateEdgeLocalDeltaTime.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [GetOffDiagonalConsistentMass](../../../../functions/containers/edge_based_data_structure/EdgeData/GetOffDiagonalConsistentMass.md)

# Called by

- [SolveSubStep](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/SolveSubStep.md)