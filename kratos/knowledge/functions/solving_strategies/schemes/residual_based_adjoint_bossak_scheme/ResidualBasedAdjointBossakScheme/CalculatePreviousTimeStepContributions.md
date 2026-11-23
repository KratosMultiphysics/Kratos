---
type: C++ Method
title: CalculatePreviousTimeStepContributions
resource: solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h#L801-L830
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
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateSystemContributions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculatePreviousTimeStepContributions( Element& rCurrentElement, LocalSystemMatrixType& rLHS_Contribution, LocalSystemVectorType& rRHS_Contribution, const ProcessInfo& rCurrentProcessInfo)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ThisThread](../../../../../functions/utilities/openmp_utils/ThisThread.md)
- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [CalculateSystemContributions](../../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateSystemContributions.md)