---
type: C++ Method
title: ReadSubModelPartEntitiesIds
resource: sources/model_part_io.cpp#L6643-L6709
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReorderedGeometryId
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadSubModelPartEntitiesIds( std::string const& rModelPartName, std::unordered_set<SizeType>& rElementsIds, std::unordered_set<SizeType>& rConditionsIds, std::unordered_set<SizeType>& rConstraintIds, std::unordered_set<SizeType>& rGeometriesIds )`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)
- [ReorderedConstraintId](../../../../functions/sources/model_part_io/ModelPartIO/ReorderedConstraintId.md)
- [ReorderedGeometryId](../../../../functions/sources/model_part_io/ModelPartIO/ReorderedGeometryId.md)