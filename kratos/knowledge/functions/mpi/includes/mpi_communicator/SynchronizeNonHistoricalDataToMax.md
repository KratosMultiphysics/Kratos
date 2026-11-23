---
type: C++ Function
title: SynchronizeNonHistoricalDataToMax
resource: mpi/includes/mpi_communicator.h#L943-L958
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/includes/mpi_communicator/TransferDistributedValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool SynchronizeNonHistoricalDataToMax(Variable<double> const& ThisVariable) override`

# Calls

- [TransferDistributedValues](../../../../functions/mpi/includes/mpi_communicator/TransferDistributedValues.md)