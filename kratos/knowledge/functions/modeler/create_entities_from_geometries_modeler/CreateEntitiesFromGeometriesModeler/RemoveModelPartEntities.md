---
type: C++ Method
title: RemoveModelPartEntities
resource: modeler/create_entities_from_geometries_modeler.cpp#L79-L108
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveElementsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateEntitiesFromGeometriesModeler::RemoveModelPartEntities(ModelPart &rModelPart)`

# Calls

- [FullName](../../../../functions/includes/model_part/FullName.md)
- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [RemoveElementsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveElementsFromAllLevels.md)
- [RemoveConditionsFromAllLevels](../../../../functions/sources/model_part/ModelPart/RemoveConditionsFromAllLevels.md)