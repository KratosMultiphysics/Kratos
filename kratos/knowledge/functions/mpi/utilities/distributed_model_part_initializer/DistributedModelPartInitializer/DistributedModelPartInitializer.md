---
type: C++ Method
title: DistributedModelPartInitializer
resource: mpi/utilities/distributed_model_part_initializer.h#L52-L59
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/tests/cpp_tests/utilities/test_distributed_model_part_initializer/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_distributed_model_part_initializer/TestDistributedModelPartInitializer/test_with_real_model
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`DistributedModelPartInitializer( ModelPart& rModelPart, const DataCommunicator& rDataComm, int SourceRank) : mrModelPart(rModelPart), mrDataComm(rDataComm), mSourceRank(SourceRank)`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/mpi/tests/cpp_tests/utilities/test_distributed_model_part_initializer/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_with_real_model](../../../../../functions/mpi/tests/test_distributed_model_part_initializer/TestDistributedModelPartInitializer/test_with_real_model.md)