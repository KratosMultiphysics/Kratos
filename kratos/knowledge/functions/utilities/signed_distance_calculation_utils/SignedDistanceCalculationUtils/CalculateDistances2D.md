---
type: C++ Method
title: CalculateDistances2D
resource: utilities/signed_distance_calculation_utils.h#L126-L336
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDistances2D( ModelPart& r_model_part, Variable<double>& rDistanceVar, const double max_distance )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [CalculateDistances](../../../../functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances.md)