---
type: Python Module
title: test_mpi_communicator
resource: mpi/tests/test_mpi_communicator.py#L1-L223
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-mpi-as-kratosmpi
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratosunittest-as-kratosunittest
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-testing-utilities
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratos-utilities
    resolved_by: tree-sitter
    confidence: exact
  - target: external/os
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GetFilePath](../../../functions/mpi/tests/test_mpi_communicator/GetFilePath.md)
- [TestMPICommunicator](../../../classes/mpi/tests/test_mpi_communicator/TestMPICommunicator.md)
- [tearDown](../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/tearDown.md)
- [_read_model_part_mpi](../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/_read_model_part_mpi.md)
- [test_assemble_variable_in_model_part](../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_assemble_variable_in_model_part.md)
- [test_assemble_variable_in_sub_model_part](../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_assemble_variable_in_sub_model_part.md)
- [testCommunicatorRankSize](../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/testCommunicatorRankSize.md)
- [testCommunicatorReduction](../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/testCommunicatorReduction.md)
- [testCommunicatorReductionSerial](../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/testCommunicatorReductionSerial.md)
- [test_GlobalNumberOf_Methods](../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_GlobalNumberOf_Methods.md)

# Imports

- `KratosMultiphysics`
- `KratosMultiphysics.mpi as KratosMPI`
- `KratosMultiphysics.KratosUnittest as KratosUnittest`
- `KratosMultiphysics.testing.utilities`
- `KratosMultiphysics.kratos_utilities`
- `os`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)