---
type: C++ Method
title: DividePropertiesBlock
resource: sources/model_part_io.cpp#L4752-L4766
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/WriteInAllFiles
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::DividePropertiesBlock(OutputFilesContainerType& OutputFiles)`

# Calls

- [WriteInAllFiles](../../../../functions/sources/model_part_io/ModelPartIO/WriteInAllFiles.md)

# Called by

- [DivideInputToPartitionsImpl](../../../../functions/sources/model_part_io/ModelPartIO/DivideInputToPartitionsImpl.md)