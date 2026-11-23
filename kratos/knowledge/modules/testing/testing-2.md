---
type: C++ Module
title: testing
resource: testing/testing.h#L1-L66
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/gtest-gtest-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/gmock-gmock-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-expect-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-data-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/testing-test-skipped-exception-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tests-test-utilities-test-suite-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tests-test-utilities-test-environment-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tests-test-utilities-test-event-listener-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GTestMain](../../classes/testing/testing/GTestMain.md)
- [InitializeTesting](../../functions/testing/testing/GTestMain/InitializeTesting.md)

# Imports

- `<gtest/gtest.h>`
- `<gmock/gmock.h>`
- `"includes/expect.h"`
- `"includes/data_communicator.h"`
- `"testing/test_skipped_exception.h"`
- `"tests/test_utilities/test_suite.h"`
- `"tests/test_utilities/test_environment.h"`
- `"tests/test_utilities/test_event_listener.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)