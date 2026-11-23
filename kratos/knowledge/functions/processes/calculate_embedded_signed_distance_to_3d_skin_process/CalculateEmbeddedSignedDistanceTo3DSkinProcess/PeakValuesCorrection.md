---
type: C++ Method
title: PeakValuesCorrection
resource: processes/calculate_embedded_signed_distance_to_3d_skin_process.h#L220-L234
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/SetMaximumAndMinimumDistanceValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PeakValuesCorrection()`

# Calls

- [SetMaximumAndMinimumDistanceValues](../../../../functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/SetMaximumAndMinimumDistanceValues.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)

# Called by

- [Execute](../../../../functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/Execute.md)