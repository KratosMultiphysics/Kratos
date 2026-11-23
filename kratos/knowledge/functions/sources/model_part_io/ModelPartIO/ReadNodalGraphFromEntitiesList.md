---
type: C++ Method
title: ReadNodalGraphFromEntitiesList
resource: sources/model_part_io.cpp#L1014-L1070
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part_io/ModelPartIO/ScanNodeBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/SkipBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlockInList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlockInList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromGeometryBlockInList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlockInList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t ModelPartIO::ReadNodalGraphFromEntitiesList( ConnectivitiesContainerType& rAuxConnectivities, std::unordered_set<SizeType>& rElementsIds, std::unordered_set<SizeType>& rConditionsIds )`

# Calls

- [ScanNodeBlock](../../../../functions/sources/model_part_io/ModelPartIO/ScanNodeBlock.md)
- [SkipBlock](../../../../functions/sources/model_part_io/ModelPartIO/SkipBlock.md)
- [FillNodalConnectivitiesFromElementBlockInList](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlockInList.md)
- [FillNodalConnectivitiesFromConditionBlockInList](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlockInList.md)
- [FillNodalConnectivitiesFromGeometryBlockInList](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromGeometryBlockInList.md)
- [FillNodalConnectivitiesFromConstraintBlockInList](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConstraintBlockInList.md)