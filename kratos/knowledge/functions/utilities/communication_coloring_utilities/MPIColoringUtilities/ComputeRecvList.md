---
type: C++ Method
title: ComputeRecvList
resource: utilities/communication_coloring_utilities.cpp#L26-L54
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/tests/cpp_tests/utilities/test_mpi_coloring_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<int> MPIColoringUtilities::ComputeRecvList( const std::vector<int>& rLocalDestinationIds, const DataCommunicator& rComm )`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/utilities/test_mpi_coloring_utilities/KRATOS_TEST_CASE_IN_SUITE.md)