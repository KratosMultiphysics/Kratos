---
type: C++ Method
title: AddVariable
resource: input_output/vtu_output.cpp#L868-L915
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__OutputTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtu_output_process/VtuOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/MakeVtuOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/WriteVtu
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtuOutput::AddVariable( const Variable<TDataType>& rVariable, const Globals::DataLocation& DataLocation)`

# Called by

- [__OutputTest](../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__OutputTest.md)
- [ExecuteInitialize](../../../../functions/python_scripts/vtu_output_process/VtuOutputProcess/ExecuteInitialize.md)
- [MakeVtuOutput](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/MakeVtuOutput.md)
- [WriteVtu](../../../../functions/tests/test_vtu_output/TestVtuOutputBase/WriteVtu.md)
- [test_PointVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition.md)
- [test_CellVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition.md)