---
type: C++ Function
title: CheckKeyIteration
resource: tests/cpp_tests/sources/test_registry.cpp#L42-L53
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/registry_item/RegistryItem/KeyConstBegin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/registry_item/RegistryItem/KeyConstEnd
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/sources/test_registry/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CheckKeyIteration( RegistryItem& rRegistryItem, const std::vector<std::string>& rRegistryItemKeys)`

# Calls

- [KeyConstBegin](../../../../../functions/sources/registry_item/RegistryItem/KeyConstBegin.md)
- [KeyConstEnd](../../../../../functions/sources/registry_item/RegistryItem/KeyConstEnd.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/sources/test_registry/KRATOS_TEST_CASE_IN_SUITE.md)