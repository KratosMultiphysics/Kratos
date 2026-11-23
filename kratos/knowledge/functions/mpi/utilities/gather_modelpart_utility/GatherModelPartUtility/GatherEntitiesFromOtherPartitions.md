---
type: C++ Method
title: GatherEntitiesFromOtherPartitions
resource: mpi/utilities/gather_modelpart_utility.cpp#L186-L220
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/tests/cpp_tests/utilities/test_gather_modelpart_utility/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GatherModelPartUtility::GatherEntitiesFromOtherPartitions( ModelPart& rModelPart, const std::map<int, std::vector<std::size_t>>& rNodesToBring, const std::map<int, std::vector<std::size_t>>& rElementsToBring, const std::map<int, std::vector<std::size_t>>& rConditionsToBring, const bool CallExecuteAfterBringingEntities, const int EchoLevel )`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/mpi/tests/cpp_tests/utilities/test_gather_modelpart_utility/KRATOS_TEST_CASE_IN_SUITE.md)