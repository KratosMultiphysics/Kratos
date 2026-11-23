---
type: C++ Method
title: ReadSubModelPartElementsAndConditionsIds
resource: sources/model_part_io.cpp#L6599-L6641
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/CheckEndBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPartIO::ReadSubModelPartElementsAndConditionsIds( std::string const& rModelPartName, std::unordered_set<SizeType>& rElementsIds, std::unordered_set<SizeType>& rConditionsIds )`

# Calls

- [CheckEndBlock](../../../../functions/sources/model_part_io/ModelPartIO/CheckEndBlock.md)