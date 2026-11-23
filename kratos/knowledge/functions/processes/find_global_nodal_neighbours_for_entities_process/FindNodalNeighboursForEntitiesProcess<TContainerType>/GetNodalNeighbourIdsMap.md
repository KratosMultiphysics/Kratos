---
type: C++ Method
title: GetNodalNeighbourIdsMap
resource: processes/find_global_nodal_neighbours_for_entities_process.cpp#L252-L309
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/find_global_nodal_neighbours_for_entities_process/GetNeighbourIds
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::unordered_map<int, std::vector<int>> FindNodalNeighboursForEntitiesProcess<TContainerType>::GetNodalNeighbourIdsMap( ModelPart::NodesContainerType& rNodes, const DataCommunicator& rDataCommunicator, const Variable<GlobalPointersVector<NodeType>>& rNodalGlobalPointerVariable)`

# Called by

- [GetNeighbourIds](../../../../functions/processes/find_global_nodal_neighbours_for_entities_process/GetNeighbourIds.md)