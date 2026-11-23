---
type: C++ Method
title: CalculateNodeNeighbourCount
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L929-L944
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetNonHistoricalVariableToZero
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Predict
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateNodeNeighbourCount(ModelPart& rModelPart)`

# Calls

- [SetNonHistoricalVariableToZero](../../../../../functions/utilities/variable_utils/SetNonHistoricalVariableToZero.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [Predict](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Predict.md)
- [InitializeSolutionStep](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/InitializeSolutionStep.md)