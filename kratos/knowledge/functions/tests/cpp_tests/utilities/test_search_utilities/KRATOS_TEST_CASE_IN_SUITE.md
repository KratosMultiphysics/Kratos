---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/utilities/test_search_utilities.cpp#L26-L41
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/utilities/test_search_utilities/GetBBoxValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/DataCommunicator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/SynchronousPointSynchronization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/SynchronousPointSynchronizationWithRadius
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(SearchUtilitiesPointIsInsideBoundingBox, KratosCoreFastSuite)`

# Calls

- [GetBBoxValue](../../../../../functions/tests/cpp_tests/utilities/test_search_utilities/GetBBoxValue.md)
- [DataCommunicator](../../../../../functions/includes/data_communicator/DataCommunicator.md)
- [SynchronousPointSynchronization](../../../../../functions/utilities/search_utilities/SynchronousPointSynchronization.md)
- [SynchronousPointSynchronizationWithRadius](../../../../../functions/utilities/search_utilities/SynchronousPointSynchronizationWithRadius.md)