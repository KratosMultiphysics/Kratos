---
type: C++ Function
title: GetNonHistoricalVariablesNames
resource: python/add_model_part_to_python.cpp#L701-L725
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/add_model_part_to_python/AddModelPartToPython
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`const std::unordered_set<std::string> GetNonHistoricalVariablesNames(ModelPart& rModelPart, TContainerType& rContainer, bool doFullSearch=false)`

# Calls

- [GetData](../../../functions/sources/gid_io/GidIOBase/GetData.md)

# Called by

- [AddModelPartToPython](../../../functions/python/add_model_part_to_python/AddModelPartToPython.md)