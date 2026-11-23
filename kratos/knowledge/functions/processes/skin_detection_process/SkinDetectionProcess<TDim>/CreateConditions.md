---
type: C++ Method
title: CreateConditions
resource: processes/skin_detection_process.cpp#L279-L315
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/RecursivelyHasProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SkinDetectionProcess<TDim>::CreateConditions( ModelPart& rMainModelPart, ModelPart& rSkinModelPart, HashMapVectorIntType& rInverseFaceMap, HashMapVectorIntIdsType& rPropertiesFaceMap, std::unordered_set<IndexType>& rNodesInTheSkin, const std::string& rConditionName) const`

# Calls

- [RecursivelyHasProperties](../../../../functions/sources/model_part/ModelPart/RecursivelyHasProperties.md)
- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)