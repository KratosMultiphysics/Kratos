---
type: C++ Method
title: ReadElementalDataBlock
resource: sources/model_part_io.cpp#L2970-L3000
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementalScalarVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementalVectorialVariableData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadInitialValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadElementalDataBlock(ElementsContainerType& rThisElements)`

# Calls

- [ReadElementalScalarVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementalScalarVariableData.md)
- [ReadElementalVectorialVariableData](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementalVectorialVariableData.md)

# Called by

- [ReadInitialValues](../../../../functions/sources/model_part_io/ModelPartIO/ReadInitialValues.md)
- [ReadModelPart](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)