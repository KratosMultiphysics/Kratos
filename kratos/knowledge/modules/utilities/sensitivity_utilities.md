---
type: C++ Module
title: sensitivity_utilities
resource: utilities/sensitivity_utilities.cpp#L1-L338
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/unordered-map
    resolved_by: tree-sitter
    confidence: exact
  - target: external/vector
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-pointer-map-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-variable-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/sensitivity-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [AddMatrixSubBlock](../../functions/utilities/sensitivity_utilities/SensitivityUtilities/AddMatrixSubBlock.md)
- [GetMatrixSubBlock](../../functions/utilities/sensitivity_utilities/SensitivityUtilities/GetMatrixSubBlock.md)
- [ComputeEntityGeometryNeighbourNodeMap](../../functions/utilities/sensitivity_utilities/SensitivityUtilities/ComputeEntityGeometryNeighbourNodeMap.md)
- [AssignEntityDerivativesToNodes](../../functions/utilities/sensitivity_utilities/SensitivityUtilities/AssignEntityDerivativesToNodes.md)

# Imports

- `<unordered_map>`
- `<vector>`
- `"includes/define.h"`
- `"includes/model_part.h"`
- `"utilities/pointer_map_communicator.h"`
- `"utilities/variable_utils.h"`
- `"utilities/parallel_utilities.h"`
- `"sensitivity_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)