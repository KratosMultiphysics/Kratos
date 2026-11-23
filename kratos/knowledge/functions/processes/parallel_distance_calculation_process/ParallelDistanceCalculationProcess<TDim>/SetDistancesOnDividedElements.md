---
type: C++ Method
title: SetDistancesOnDividedElements
resource: processes/parallel_distance_calculation_process.cpp#L311-L356
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/IsDivided
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/PreserveDistancesOnDividedElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/CalculateExactDistancesOnDividedElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelDistanceCalculationProcess<TDim>::SetDistancesOnDividedElements()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [IsDivided](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/IsDivided.md)
- [PreserveDistancesOnDividedElements](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/PreserveDistancesOnDividedElements.md)
- [CalculateExactDistancesOnDividedElements](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/CalculateExactDistancesOnDividedElements.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)

# Called by

- [Execute](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/Execute.md)