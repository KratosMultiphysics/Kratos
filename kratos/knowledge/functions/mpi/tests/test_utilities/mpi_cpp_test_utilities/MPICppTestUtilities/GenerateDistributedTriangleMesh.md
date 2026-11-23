---
type: C++ Method
title: GenerateDistributedTriangleMesh
resource: mpi/tests/test_utilities/mpi_cpp_test_utilities.cpp#L185-L271
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/cpp_tests/utilities/test_parallel_fill_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPICppTestUtilities::GenerateDistributedTriangleMesh( ModelPart& rModelPart, const DataCommunicator& rDataCommunicator )`

# Calls

- [CreateNewProperties](../../../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../../functions/mpi/tests/cpp_tests/utilities/test_parallel_fill_communicator/KRATOS_TEST_CASE_IN_SUITE.md)