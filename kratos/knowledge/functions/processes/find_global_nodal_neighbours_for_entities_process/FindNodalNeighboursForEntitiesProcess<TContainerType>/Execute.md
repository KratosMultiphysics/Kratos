---
type: C++ Method
title: Execute
resource: processes/find_global_nodal_neighbours_for_entities_process.cpp#L76-L234
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/communication_coloring_utilities/MPIColoringUtilities/ComputeCommunicationScheduling
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindNodalNeighboursForEntitiesProcess<TContainerType>::Execute()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ComputeCommunicationScheduling](../../../../functions/utilities/communication_coloring_utilities/MPIColoringUtilities/ComputeCommunicationScheduling.md)
- [RetrieveGlobalIndexedPointersMap](../../../../functions/utilities/global_pointer_utilities/RetrieveGlobalIndexedPointersMap.md)