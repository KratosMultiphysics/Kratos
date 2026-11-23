---
type: C++ Method
title: CalculateDistanceByBruteForce
resource: processes/calculate_distance_to_path_process.cpp#L202-L229
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/FastMinimalDistanceOnLine
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/GeometricalProjectionUtilities/FastMinimalDistanceOnLineWithRadius
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDistanceToPathProcess<THistorical>::CalculateDistanceByBruteForce( ModelPart& rModelPart, std::vector<Geometry<NodeType>::Pointer>& rVectorSegments )`

# Calls

- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [FastMinimalDistanceOnLine](../../../../functions/utilities/geometrical_projection_utilities/FastMinimalDistanceOnLine.md)
- [FastMinimalDistanceOnLineWithRadius](../../../../functions/utilities/geometrical_projection_utilities/GeometricalProjectionUtilities/FastMinimalDistanceOnLineWithRadius.md)

# Called by

- [Execute](../../../../functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/Execute.md)