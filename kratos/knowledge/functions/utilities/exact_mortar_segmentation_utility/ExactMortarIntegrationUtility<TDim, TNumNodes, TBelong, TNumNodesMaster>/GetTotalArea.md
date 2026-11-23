---
type: C++ Method
title: GetTotalArea
resource: utilities/exact_mortar_segmentation_utility.cpp#L1047-L1067
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/GetExactAreaIntegration
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::GetTotalArea( const GeometryType& rOriginalSlaveGeometry, ConditionArrayListType& rConditionsPointsSlave, double& rArea )`

# Called by

- [GetExactAreaIntegration](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/GetExactAreaIntegration.md)