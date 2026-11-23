---
type: C++ Method
title: GetExactIntegration
resource: utilities/exact_mortar_segmentation_utility.cpp#L979-L1020
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::GetExactIntegration( const GeometryType& rOriginalSlaveGeometry, const array_1d<double, 3>& rSlaveNormal, const GeometryType& rOriginalMasterGeometry, const array_1d<double, 3>& rMasterNormal, IntegrationPointsType& rIntegrationPointsSlave )`

# Calls

- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)