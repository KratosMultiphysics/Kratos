---
type: C++ Method
title: DivideFlagVariableData
resource: sources/model_part_io.cpp#L5220-L5260
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/DivideNodalDataBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::DivideFlagVariableData(OutputFilesContainerType& OutputFiles, PartitionIndicesContainerType const& NodesAllPartitions)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [DivideNodalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideNodalDataBlock.md)