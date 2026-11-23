---
type: C++ Method
title: CalculateAe
resource: processes/simple_mortar_mapper_process.cpp#L452-L503
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/mortar_utilities/LengthCheck
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAeComponents
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline BoundedMatrix<double, TNumNodes, TNumNodes> SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::CalculateAe( GeometryType& rSlaveGeometry, MortarKinematicVariablesType& rThisKinematicVariables, std::vector<array_1d<PointType,TDim>>& rGeometricalObjectsPointSlave, const IntegrationMethod& rThisIntegrationMethod )`

# Calls

- [LengthCheck](../../../../functions/utilities/mortar_utilities/LengthCheck.md)
- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)
- [CalculateAeComponents](../../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAeComponents.md)