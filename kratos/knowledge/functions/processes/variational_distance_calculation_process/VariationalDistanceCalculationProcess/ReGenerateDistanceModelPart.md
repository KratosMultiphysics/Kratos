---
type: C++ Method
title: ReGenerateDistanceModelPart
resource: processes/variational_distance_calculation_process.h#L557-L597
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/model/Model/HasModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/DeleteModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/VariationalDistanceCalculationProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void ReGenerateDistanceModelPart(ModelPart& rBaseModelPart)`

# Calls

- [HasModelPart](../../../../functions/containers/model/Model/HasModelPart.md)
- [DeleteModelPart](../../../../functions/containers/model/Model/DeleteModelPart.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [VariationalDistanceCalculationProcess](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/VariationalDistanceCalculationProcess.md)
- [Execute](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute.md)