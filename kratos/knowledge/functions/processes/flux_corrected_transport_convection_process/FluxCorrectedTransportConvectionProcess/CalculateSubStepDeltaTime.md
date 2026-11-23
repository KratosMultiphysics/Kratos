---
type: C++ Method
title: CalculateSubStepDeltaTime
resource: processes/flux_corrected_transport_convection_process.h#L314-L357
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateEdgeLocalDeltaTime
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double CalculateSubStepDeltaTime()`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [CalculateEdgeLocalDeltaTime](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateEdgeLocalDeltaTime.md)

# Called by

- [Execute](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/Execute.md)