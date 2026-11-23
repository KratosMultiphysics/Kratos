---
type: C++ Method
title: SetMaximumAndMinimumDistanceValues
resource: processes/calculate_embedded_signed_distance_to_3d_skin_process.h#L236-L293
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/PeakValuesCorrection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetMaximumAndMinimumDistanceValues(double& max_distance, double& min_distance)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)

# Called by

- [PeakValuesCorrection](../../../../functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/PeakValuesCorrection.md)