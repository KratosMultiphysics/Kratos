---
type: C++ Function
title: SynchronizeNonHistoricalDataToMin
resource: mpi/includes/mpi_communicator.h#L1011-L1026
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/includes/mpi_communicator/TransferDistributedValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool SynchronizeNonHistoricalDataToMin(Variable<double> const& ThisVariable) override`

# Calls

- [TransferDistributedValues](../../../../functions/mpi/includes/mpi_communicator/TransferDistributedValues.md)