---
type: C++ Function
title: SynchronousPointSynchronization
resource: utilities/search_utilities.h#L311-L328
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/search_utilities/CalculateNumberOfPoints
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/cpp_tests/utilities/test_search_utilities/KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_search_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void SynchronousPointSynchronization( TPointIteratorType itPointBegin, TPointIteratorType itPointEnd, std::vector<double>& rAllPointsCoordinates, std::vector<SignedIndexType>& rAllPointsIds, const DataCommunicator& rDataCommunicator )`

# Calls

- [CalculateNumberOfPoints](../../../functions/utilities/search_utilities/CalculateNumberOfPoints.md)

# Called by

- [KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/utilities/test_search_utilities/KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_search_utilities/KRATOS_TEST_CASE_IN_SUITE.md)