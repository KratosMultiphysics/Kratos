---
type: C++ Method
title: DivideElementsBlock
resource: sources/model_part_io.cpp#L4895-L4970
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteInAllFiles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::DivideElementsBlock(OutputFilesContainerType& OutputFiles, PartitionIndicesContainerType const& ElementsAllPartitions)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [WriteInAllFiles](../../../../functions/sources/model_part_io/ModelPartIO/WriteInAllFiles.md)
- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [DivideInputToPartitionsImpl](../../../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl.md)