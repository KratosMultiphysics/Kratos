---
type: C++ Method
title: CalculateDistances3D
resource: utilities/signed_distance_calculation_utils.h#L338-L456
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateTetrahedraDistances
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDistances3D( ModelPart& r_model_part, Variable<double>& rDistanceVar, const double max_distance )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CalculateTetrahedraDistances](../../../../functions/utilities/geometry_utilities/CalculateTetrahedraDistances.md)

# Called by

- [CalculateDistances](../../../../functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances.md)