---
type: C++ Method
title: MPIIsFinalizedDetail
resource: sources/parallel_environment.cpp#L359-L362
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/sources/mpi_manager/MPIManager/IsFinalized
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/parallel_environment/ParallelEnvironment/MPIIsFinalized
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ParallelEnvironment::MPIIsFinalizedDetail() const`

# Calls

- [IsFinalized](../../../../functions/mpi/sources/mpi_manager/MPIManager/IsFinalized.md)

# Called by

- [MPIIsFinalized](../../../../functions/sources/parallel_environment/ParallelEnvironment/MPIIsFinalized.md)