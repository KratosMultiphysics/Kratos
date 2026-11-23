---
type: C++ Method
title: Execute
resource: processes/variational_distance_calculation_process.h#L284-L421
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/ReGenerateDistanceModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Free
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateExactDistancesToPlane
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateTetrahedraDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateTriangleDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/SynchronizeFixity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/SynchronizeDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Execute() override`

# Calls

- [ReGenerateDistanceModelPart](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/ReGenerateDistanceModelPart.md)
- [Free](../../../../functions/includes/node/Free.md)
- [Fix](../../../../functions/includes/node/Fix.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CalculateExactDistancesToPlane](../../../../functions/utilities/geometry_utilities/CalculateExactDistancesToPlane.md)
- [CalculateTetrahedraDistances](../../../../functions/utilities/geometry_utilities/CalculateTetrahedraDistances.md)
- [CalculateTriangleDistances](../../../../functions/utilities/geometry_utilities/CalculateTriangleDistances.md)
- [SynchronizeFixity](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/SynchronizeFixity.md)
- [SynchronizeDistance](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/SynchronizeDistance.md)
- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)