---
type: C++ Method
title: VtuOutput
resource: input_output/vtu_output.cpp#L753-L785
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtu_output/AddUnstructuredGridData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__OutputTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtu_output_process/VtuOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/WriteVtu
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/test_GetOutputContainerList
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

`VtuOutput::VtuOutput( ModelPart& rModelPart, const bool IsInitialConfiguration, const WriterFormat OutputFormat, const IndexType Precision, const bool OutputSubModelParts, const bool WriteIds, const IndexType EchoLevel) : mIsSeriesFileHeaderWritten(false), mrModelPart(rModelPart), mConfiguration(IsInitialConfiguration ? Globals::Configuration::Initial : Globals::Configuration::Current), mEchoLevel(EchoLevel), mOutputFormat(OutputFormat), mPrecision(Precision)`

# Calls

- [AddUnstructuredGridData](../../../../functions/input_output/vtu_output/AddUnstructuredGridData.md)
- [FullName](../../../../functions/includes/model_part/FullName.md)

# Called by

- [__OutputTest](../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__OutputTest.md)
- [ExecuteInitialize](../../../../functions/python_scripts/vtu_output_process/VtuOutputProcess/ExecuteInitialize.md)
- [WriteVtu](../../../../functions/tests/test_vtu_output/TestVtuOutputBase/WriteVtu.md)
- [test_GetOutputContainerList](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_GetOutputContainerList.md)
- [test_PointVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_PointVariableAddition.md)
- [test_CellVariableAddition](../../../../functions/tests/test_vtu_output/TestVtuOutput/test_CellVariableAddition.md)