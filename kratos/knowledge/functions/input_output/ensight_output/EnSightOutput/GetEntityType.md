---
type: C++ Method
title: GetEntityType
resource: input_output/ensight_output.cpp#L259-L282
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/GlobalNumberOfElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfConditions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`EnSightOutput::EntityType EnSightOutput::GetEntityType(const ModelPart& rModelPart) const`

# Calls

- [GlobalNumberOfElements](../../../../functions/sources/communicator/Communicator/GlobalNumberOfElements.md)
- [GlobalNumberOfConditions](../../../../functions/sources/communicator/Communicator/GlobalNumberOfConditions.md)