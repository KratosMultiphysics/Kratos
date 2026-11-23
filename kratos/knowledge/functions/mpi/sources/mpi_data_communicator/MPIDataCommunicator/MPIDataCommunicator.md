---
type: C++ Method
title: MPIDataCommunicator
resource: mpi/sources/mpi_data_communicator.cpp#L293-L301
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/MPIIsInitialized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/parallel_environment/ParallelEnvironment/SetUpMPIEnvironment
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`MPIDataCommunicator::MPIDataCommunicator(MPI_Comm MPIComm): DataCommunicator(), mComm(MPIComm)`

# Calls

- [MPIIsInitialized](../../../../../functions/sources/parallel_environment/ParallelEnvironment/MPIIsInitialized.md)
- [SetUpMPIEnvironment](../../../../../functions/sources/parallel_environment/ParallelEnvironment/SetUpMPIEnvironment.md)