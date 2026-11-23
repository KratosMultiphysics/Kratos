---
type: C++ Function
title: UpdateEntityTimeSchemeContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L1006-L1047
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/openmp_utils/ThisThread
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateTimeSchemeAdjoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void UpdateEntityTimeSchemeContributions( TEntityContainerType& rEntityContainer, const ProcessInfo& rProcessInfo)`

# Calls

- [ThisThread](../../../../functions/utilities/openmp_utils/ThisThread.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [UpdateTimeSchemeAdjoints](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/UpdateTimeSchemeAdjoints.md)