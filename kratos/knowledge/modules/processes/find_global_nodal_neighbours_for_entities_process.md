---
type: C++ Module
title: find_global_nodal_neighbours_for_entities_process
resource: processes/find_global_nodal_neighbours_for_entities_process.cpp#L1-L326
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/string
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-element-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-global-pointer-variables-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-node-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-communication-coloring-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-global-pointer-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-pointer-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-variable-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/find-global-nodal-neighbours-for-entities-process-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [AddHangingNodeIds](../../functions/processes/find_global_nodal_neighbours_for_entities_process/FindNodalNeighboursForEntitiesProcess<ModelPart/ElementsContainerType>/AddHangingNodeIds.md)
- [AddHangingNodeIds](../../functions/processes/find_global_nodal_neighbours_for_entities_process/FindNodalNeighboursForEntitiesProcess<ModelPart/ConditionsContainerType>/AddHangingNodeIds.md)
- [Execute](../../functions/processes/find_global_nodal_neighbours_for_entities_process/FindNodalNeighboursForEntitiesProcess<TContainerType>/Execute.md)
- [ClearNeighbours](../../functions/processes/find_global_nodal_neighbours_for_entities_process/FindNodalNeighboursForEntitiesProcess<TContainerType>/ClearNeighbours.md)
- [GetNodalNeighbourIdsMap](../../functions/processes/find_global_nodal_neighbours_for_entities_process/FindNodalNeighboursForEntitiesProcess<TContainerType>/GetNodalNeighbourIdsMap.md)
- [GlobalPointerAdder](../../classes/processes/find_global_nodal_neighbours_for_entities_process/GlobalPointerAdder.md)
- [GetValue](../../functions/processes/find_global_nodal_neighbours_for_entities_process/GlobalPointerAdder/GetValue.md)
- [LocalReduce](../../functions/processes/find_global_nodal_neighbours_for_entities_process/GlobalPointerAdder/LocalReduce.md)
- [ThreadSafeReduce](../../functions/processes/find_global_nodal_neighbours_for_entities_process/GlobalPointerAdder/ThreadSafeReduce.md)
- [AddUnique](../../functions/processes/find_global_nodal_neighbours_for_entities_process/FindNodalNeighboursForEntitiesProcess<TContainerType>/AddUnique.md)

# Imports

- `<iostream>`
- `<string>`
- `"includes/define.h"`
- `"includes/element.h"`
- `"includes/global_pointer_variables.h"`
- `"includes/model_part.h"`
- `"includes/node.h"`
- `"processes/process.h"`
- `"utilities/communication_coloring_utilities.h"`
- `"utilities/global_pointer_utilities.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/pointer_communicator.h"`
- `"utilities/variable_utils.h"`
- `"find_global_nodal_neighbours_for_entities_process.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)