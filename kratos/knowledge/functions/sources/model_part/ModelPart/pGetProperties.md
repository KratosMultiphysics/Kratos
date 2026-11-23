---
type: C++ Method
title: pGetProperties
resource: sources/model_part.cpp#L694-L714
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/properties/HasSubProperties
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPart::PropertiesType::Pointer ModelPart::pGetProperties( IndexType PropertiesId, IndexType MeshIndex )`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [HasSubProperties](../../../../functions/includes/properties/HasSubProperties.md)