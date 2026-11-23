---
type: C++ Method
title: TotalProcesses
resource: sources/communicator.cpp#L96-L99
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtu_output/VtuOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/testCommunicatorRankSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/SynchronizeDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/SynchronizeFixity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/json_output_process/JsonOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/PrintData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int Communicator::TotalProcesses() const`

# Called by

- [PrintOutput](../../../../functions/input_output/vtu_output/VtuOutput/PrintOutput.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE.md)
- [testCommunicatorRankSize](../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/testCommunicatorRankSize.md)
- [SynchronizeDistance](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/SynchronizeDistance.md)
- [SynchronizeFixity](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/SynchronizeFixity.md)
- [ExecuteInitialize](../../../../functions/python_scripts/json_output_process/JsonOutputProcess/ExecuteInitialize.md)
- [__init__](../../../../functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__init__.md)
- [PrintData](../../../../functions/sources/model_part/ModelPart/PrintData.md)