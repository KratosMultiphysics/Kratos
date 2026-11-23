---
type: C++ Method
title: SetUpMPIEnvironment
resource: sources/parallel_environment.cpp#L56-L60
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/parallel_environment/ParallelEnvironment/SetUpMPIEnvironmentDetail
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/MPIDataCommunicator
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ParallelEnvironment::SetUpMPIEnvironment(EnvironmentManager::Pointer pEnvironmentManager)`

# Calls

- [SetUpMPIEnvironmentDetail](../../../../functions/sources/parallel_environment/ParallelEnvironment/SetUpMPIEnvironmentDetail.md)

# Called by

- [MPIDataCommunicator](../../../../functions/mpi/sources/mpi_data_communicator/MPIDataCommunicator/MPIDataCommunicator.md)