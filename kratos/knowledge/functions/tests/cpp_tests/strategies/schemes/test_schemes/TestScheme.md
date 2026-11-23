---
type: C++ Function
title: TestScheme
resource: tests/cpp_tests/strategies/schemes/test_schemes.cpp#L113-L385
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/strategies/schemes/test_schemes/BasicTestSchemeDisplacement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FreeDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FixDof
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/strategies/schemes/test_schemes/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void TestScheme( SchemeType::Pointer pScheme, const double DeltaTime = 1.0e-4, const std::string TestType = "DISPLACEMENT", const bool TestPredict = false, const bool TestUpdatePredict = false )`

# Calls

- [BasicTestSchemeDisplacement](../../../../../../functions/tests/cpp_tests/strategies/schemes/test_schemes/BasicTestSchemeDisplacement.md)
- [FreeDof](../../../../../../functions/includes/dof/Dof/FreeDof.md)
- [FixDof](../../../../../../functions/includes/dof/Dof/FixDof.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../../functions/tests/cpp_tests/strategies/schemes/test_schemes/KRATOS_TEST_CASE_IN_SUITE.md)