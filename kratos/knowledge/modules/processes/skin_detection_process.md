---
type: C++ Module
title: skin_detection_process
resource: processes/skin_detection_process.cpp#L1-L581
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/processes-skin-detection-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-variable-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-key-hash-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-reduction-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [SkinDetectionProcess](../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/SkinDetectionProcess.md)
- [Execute](../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/Execute.md)
- [GenerateFaceMaps](../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/GenerateFaceMaps.md)
- [FillAuxiliaryModelPart](../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/FillAuxiliaryModelPart.md)
- [CreateConditions](../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/CreateConditions.md)
- [SetUpAdditionalSubModelParts](../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/SetUpAdditionalSubModelParts.md)
- [GenerateSetNodeIdsInterface](../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/GenerateSetNodeIdsInterface.md)
- [FilterMPIInterfaceNodes](../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/FilterMPIInterfaceNodes.md)
- [GetDefaultParameters](../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/GetDefaultParameters.md)
- [GetSettings](../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/GetSettings.md)
- [SkinDetectionProcess](../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/SkinDetectionProcess-2.md)

# Imports

- `"processes/skin_detection_process.h"`
- `"utilities/variable_utils.h"`
- `"includes/key_hash.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/reduction_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)