---
type: C++ Method
title: ReplaceTensorAdaptor
resource: input_output/vtu_output.cpp#L1020-L1072
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtu_output/FindUnstructuredGridData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtu_output/VtuOutput/EmplaceTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_mpi_vtu_output/AddTensorAdaptor
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtuOutput::ReplaceTensorAdaptor( const std::string& rTensorAdaptorName, TTensorAdaptorPointerType pTensorAdaptor)`

# Calls

- [FindUnstructuredGridData](../../../../functions/input_output/vtu_output/FindUnstructuredGridData.md)

# Called by

- [EmplaceTensorAdaptor](../../../../functions/input_output/vtu_output/VtuOutput/EmplaceTensorAdaptor.md)
- [AddTensorAdaptor](../../../../functions/mpi/future/tests/test_mpi_vtu_output/AddTensorAdaptor.md)