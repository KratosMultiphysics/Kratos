---
type: C++ Method
title: UpdateSensitivities
resource: utilities/sensitivity_builder.cpp#L418-L460
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/sensitivity_builder/SensitivityBuilder/ClearSensitivities
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/RunAdjointSensitivityTest
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SensitivityBuilder::UpdateSensitivities()`

# Calls

- [ClearSensitivities](../../../../functions/utilities/sensitivity_builder/SensitivityBuilder/ClearSensitivities.md)

# Called by

- [RunAdjointSensitivityTest](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/RunAdjointSensitivityTest.md)