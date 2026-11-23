---
type: C++ Function
title: BasicTestSchemeDisplacement
resource: tests/cpp_tests/strategies/schemes/test_schemes.cpp#L45-L108
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/strategies/schemes/test_schemes/TestScheme
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static DofsArrayType BasicTestSchemeDisplacement( ModelPart& rModelPart, SchemeType::Pointer pScheme, std::vector< Dof<double>::Pointer >& rDoF, const double DeltaTime = 1.0e-4 )`

# Called by

- [TestScheme](../../../../../../functions/tests/cpp_tests/strategies/schemes/test_schemes/TestScheme.md)