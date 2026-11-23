---
type: C++ Method
title: MPIIsInitialized
resource: sources/parallel_environment.cpp#L141-L145
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/MPIIsInitializedDetail
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/MPIDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/parallel_environment/ParallelEnvironment/SetUpMPIEnvironmentDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ParallelEnvironment::MPIIsInitialized()`

# Calls

- [MPIIsInitializedDetail](../../../../functions/sources/parallel_environment/ParallelEnvironment/MPIIsInitializedDetail.md)

# Called by

- [MPIDataCommunicator](../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/MPIDataCommunicator.md)
- [SetUpMPIEnvironmentDetail](../../../../functions/sources/parallel_environment/ParallelEnvironment/SetUpMPIEnvironmentDetail.md)