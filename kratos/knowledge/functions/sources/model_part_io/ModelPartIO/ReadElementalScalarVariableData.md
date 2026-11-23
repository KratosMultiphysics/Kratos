---
type: C++ Method
title: ReadElementalScalarVariableData
resource: sources/model_part_io.cpp#L3003-L3032
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementalDataBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadElementalScalarVariableData(ElementsContainerType& rThisElements, const TVariableType& rVariable)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)

# Called by

- [ReadElementalDataBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementalDataBlock.md)