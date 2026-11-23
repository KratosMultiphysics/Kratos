---
type: C++ Function
title: AssembleNonHistoricalData
resource: mpi/includes/mpi_communicator.h#L1097-L1102
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/includes/mpi_communicator/AssembleFixedSizeValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/AssembleDynamicVectorValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/includes/mpi_communicator/AssembleDynamicMatrixValues
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool AssembleNonHistoricalData(Variable<int> const& ThisVariable) override`

# Calls

- [AssembleFixedSizeValues](../../../../functions/mpi/includes/mpi_communicator/AssembleFixedSizeValues.md)
- [AssembleDynamicVectorValues](../../../../functions/mpi/includes/mpi_communicator/AssembleDynamicVectorValues.md)
- [AssembleDynamicMatrixValues](../../../../functions/mpi/includes/mpi_communicator/AssembleDynamicMatrixValues.md)