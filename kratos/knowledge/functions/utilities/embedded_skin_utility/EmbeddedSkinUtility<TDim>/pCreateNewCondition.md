---
type: C++ Method
title: pCreateNewCondition
resource: utilities/embedded_skin_utility.cpp#L242-L250
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/pCreateNewConditionGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/ComputeElementSkin
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Condition::Pointer EmbeddedSkinUtility<TDim>::pCreateNewCondition( const GeometryData::KratosGeometryType &rOriginGeometryType, const Condition::NodesArrayType &rNewNodesArray, const unsigned int &rConditionId, const Properties::Pointer pConditionProperties)`

# Calls

- [pCreateNewConditionGeometry](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/pCreateNewConditionGeometry.md)

# Called by

- [ComputeElementSkin](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/ComputeElementSkin.md)