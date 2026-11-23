---
type: C++ Function
title: ModelPartForMPICommunicatorTests
resource: mpi/tests/cpp_tests/sources/test_mpi_communicator.cpp#L31-L71
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartForMPICommunicatorTests(ModelPart& rModelPart, const DataCommunicator& rComm)`

# Calls

- [CreateNewProperties](../../../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../../functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE.md)