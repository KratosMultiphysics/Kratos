---
type: C++ Function
title: RunAdjointSensitivityTest
resource: tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme.cpp#L988-L1063
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/InitializePrimalModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/InitializeAdjointModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_builder/SensitivityBuilder/UpdateSensitivities
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double RunAdjointSensitivityTest( Model& rModel, const bool IsWithConditions, const double Perturbation)`

# Calls

- [InitializePrimalModelPart](../../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/InitializePrimalModelPart.md)
- [InitializeAdjointModelPart](../../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/InitializeAdjointModelPart.md)
- [UpdateSensitivities](../../../../../../functions/utilities/sensitivity_builder/SensitivityBuilder/UpdateSensitivities.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/KRATOS_TEST_CASE_IN_SUITE.md)