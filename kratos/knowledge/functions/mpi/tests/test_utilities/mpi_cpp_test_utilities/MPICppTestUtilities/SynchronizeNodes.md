---
type: C++ Method
title: SynchronizeNodes
resource: mpi/tests/test_utilities/mpi_cpp_test_utilities.cpp#L369-L380
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/global_pointer_utilities/GlobalRetrieveGlobalPointers
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/cpp_tests/utilities/test_gather_modelpart_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GlobalPointerCommunicator<Node>::Pointer MPICppTestUtilities::SynchronizeNodes( ModelPart& rModelPart, GlobalPointersVector<Node>& rGlobalPointers, const DataCommunicator& rDataCommunicator )`

# Calls

- [GlobalRetrieveGlobalPointers](../../../../../../functions/utilities/global_pointer_utilities/GlobalRetrieveGlobalPointers.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../../functions/mpi/tests/cpp_tests/utilities/test_gather_modelpart_utility/KRATOS_TEST_CASE_IN_SUITE.md)