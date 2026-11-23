---
type: C++ Method
title: RegisterDataCommunicator
resource: sources/parallel_environment.cpp#L114-L121
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicatorDetail
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/python/kratos_mpi_python/InitializeMPIParallelRun
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_utilities/mpi_test_environment/KratosMpiTestEnv/SetUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_utilities/mpi_test_environment/KratosMpiTestEnv/TearDown
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelEnvironment::RegisterDataCommunicator( const std::string& Name, DataCommunicator::UniquePointer pPrototype, const bool Default)`

# Calls

- [RegisterDataCommunicatorDetail](../../../../functions/sources/parallel_environment/ParallelEnvironment/RegisterDataCommunicatorDetail.md)

# Called by

- [InitializeMPIParallelRun](../../../../functions/mpi/python/kratos_mpi_python/InitializeMPIParallelRun.md)
- [SetUp](../../../../functions/mpi/tests/test_utilities/mpi_test_environment/KratosMpiTestEnv/SetUp.md)
- [TearDown](../../../../functions/mpi/tests/test_utilities/mpi_test_environment/KratosMpiTestEnv/TearDown.md)