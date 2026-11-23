---
type: C++ Module
title: distributed_model_part_initializer
resource: mpi/utilities/distributed_model_part_initializer.cpp#L1-L101
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/utilities-string-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/parallel-fill-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/model-part-communicator-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/distributed-model-part-initializer-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GetSubModelPartHierarchy](../../../functions/mpi/utilities/distributed_model_part_initializer/GetSubModelPartHierarchy.md)
- [RecursiveCreateModelParts](../../../functions/mpi/utilities/distributed_model_part_initializer/RecursiveCreateModelParts.md)
- [CreateSubModelPartHierarchy](../../../functions/mpi/utilities/distributed_model_part_initializer/CreateSubModelPartHierarchy.md)
- [CopySubModelPartStructure](../../../functions/mpi/utilities/distributed_model_part_initializer/DistributedModelPartInitializer/CopySubModelPartStructure.md)
- [Execute](../../../functions/mpi/utilities/distributed_model_part_initializer/DistributedModelPartInitializer/Execute.md)

# Imports

- `"utilities/string_utilities.h"`
- `"parallel_fill_communicator.h"`
- `"model_part_communicator_utilities.h"`
- `"distributed_model_part_initializer.h"`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)