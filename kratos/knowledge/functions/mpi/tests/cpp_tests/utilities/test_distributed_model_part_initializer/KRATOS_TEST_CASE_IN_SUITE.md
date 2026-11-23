---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: mpi/tests/cpp_tests/utilities/test_distributed_model_part_initializer.cpp#L25-L37
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/utilities/distributed_model_part_initializer/DistributedModelPartInitializer/DistributedModelPartInitializer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(DistributedModelPartInitializer_NoSubModelParts_Empty, KratosMPICoreFastSuite)`

# Calls

- [DistributedModelPartInitializer](../../../../../../functions/mpi/utilities/distributed_model_part_initializer/DistributedModelPartInitializer/DistributedModelPartInitializer.md)
- [NumberOfSubModelParts](../../../../../../functions/includes/model_part/NumberOfSubModelParts.md)
- [GlobalNumberOfNodes](../../../../../../functions/sources/communicator/Communicator/GlobalNumberOfNodes.md)
- [HasSubModelPart](../../../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)