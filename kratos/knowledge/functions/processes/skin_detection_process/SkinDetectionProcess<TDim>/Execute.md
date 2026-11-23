---
type: C++ Method
title: Execute
resource: processes/skin_detection_process.cpp#L40-L62
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/GenerateSetNodeIdsInterface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/GenerateFaceMaps
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/FilterMPIInterfaceNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/FillAuxiliaryModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/SetUpAdditionalSubModelParts
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SkinDetectionProcess<TDim>::Execute()`

# Calls

- [GenerateSetNodeIdsInterface](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/GenerateSetNodeIdsInterface.md)
- [GenerateFaceMaps](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/GenerateFaceMaps.md)
- [FilterMPIInterfaceNodes](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/FilterMPIInterfaceNodes.md)
- [FillAuxiliaryModelPart](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/FillAuxiliaryModelPart.md)
- [SetUpAdditionalSubModelParts](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/SetUpAdditionalSubModelParts.md)