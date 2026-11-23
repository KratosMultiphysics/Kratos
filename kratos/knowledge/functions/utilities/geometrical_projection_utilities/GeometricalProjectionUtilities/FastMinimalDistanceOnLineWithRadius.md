---
type: C++ Method
title: FastMinimalDistanceOnLineWithRadius
resource: utilities/geometrical_projection_utilities.cpp#L24-L104
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/intersection_utilities/ComputeShortestLineBetweenTwoLines
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistanceByBruteForce
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GeometricalProjectionUtilities::DistanceComputed GeometricalProjectionUtilities::FastMinimalDistanceOnLineWithRadius( double& rDistance, const Geometry<Node>& rSegment, const Point& rPoint, const double Radius, const double Tolerance )`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [ComputeShortestLineBetweenTwoLines](../../../../functions/utilities/intersection_utilities/ComputeShortestLineBetweenTwoLines.md)

# Called by

- [CalculateDistance](../../../../functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistance.md)
- [CalculateDistanceByBruteForce](../../../../functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistanceByBruteForce.md)