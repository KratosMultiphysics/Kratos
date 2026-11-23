---
type: C++ Function
title: SynchronousPointSynchronizationWithRadius
resource: utilities/search_utilities.h#L412-L437
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/search_utilities/CalculateNumberOfPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/SynchronizeRadius
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

`static std::vector<double> SynchronousPointSynchronizationWithRadius( TPointIteratorType itPointBegin, TPointIteratorType itPointEnd, std::vector<double>& rAllPointsCoordinates, std::vector<SignedIndexType>& rAllPointsIds, const std::vector<double>& rRadius, const DataCommunicator& rDataCommunicator )`

# Calls

- [CalculateNumberOfPoints](../../../functions/utilities/search_utilities/CalculateNumberOfPoints.md)
- [SynchronizeRadius](../../../functions/utilities/search_utilities/SynchronizeRadius.md)

# Called by

- [KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/utilities/test_search_utilities/KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_search_utilities/KRATOS_TEST_CASE_IN_SUITE.md)