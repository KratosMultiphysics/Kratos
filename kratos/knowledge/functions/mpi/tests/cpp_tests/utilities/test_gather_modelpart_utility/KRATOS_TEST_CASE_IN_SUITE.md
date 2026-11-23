---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: mpi/tests/cpp_tests/utilities/test_gather_modelpart_utility.cpp#L27-L100
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/tests/test_utilities/mpi_cpp_test_utilities/MPICppTestUtilities/GenerateDistributedBarStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_utilities/mpi_cpp_test_utilities/MPICppTestUtilities/SynchronizeNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/gather_modelpart_utility/GatherModelPartUtility/GatherEntitiesFromOtherPartitions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(GatherModelPartUtilityGatherEntitiesFromOtherPartitions, KratosMPICoreFastSuite)`

# Calls

- [GenerateDistributedBarStructure](../../../../../../functions/mpi/tests/test_utilities/mpi_cpp_test_utilities/MPICppTestUtilities/GenerateDistributedBarStructure.md)
- [SynchronizeNodes](../../../../../../functions/mpi/tests/test_utilities/mpi_cpp_test_utilities/MPICppTestUtilities/SynchronizeNodes.md)
- [GatherEntitiesFromOtherPartitions](../../../../../../functions/mpi/utilities/gather_modelpart_utility/GatherModelPartUtility/GatherEntitiesFromOtherPartitions.md)