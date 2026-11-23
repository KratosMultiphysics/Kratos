---
type: C++ Method
title: ReadPropertiesBlock
resource: sources/model_part_io.cpp#L2066-L2169
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstitutiveLawValue
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/model_part_io/ModelPartIO/ReadProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadPropertiesBlock(PropertiesContainerType& rThisProperties)`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReadConstitutiveLawValue](../../../../functions/sources/model_part_io/ModelPartIO/ReadConstitutiveLawValue.md)

# Called by

- [ReadProperties](../../../../functions/sources/model_part_io/ModelPartIO/ReadProperties.md)
- [ReadModelPart](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)